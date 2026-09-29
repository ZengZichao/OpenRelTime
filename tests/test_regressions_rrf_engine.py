"""Regression gates for the core RRF engine.

Each test pins one behaviour the engine must keep satisfying:

* ``mean="arithmetic"`` folds child-subtree lengths with the *arithmetic* mean,
  so the convention is never mixed and ``b = r * t`` holds on trees deeper than
  two levels;
* the rate-ratio guard's registered default is actually exercised: these trees
  do trigger it, and both the replaced node ages and *which* ancestor supplied
  the replacement are pinned, not merely "some ancestor";
* ``rrf_rates`` reports honestly that it applies no guard, because a rate-only
  pass has no ages to guard;
* ``RrfEngine._composite`` is the only folding path, and every folding site
  routes through it;
* the guard defaults come from ``openreltime._constants``, the single source of
  registered defaults (see docs/parameters.md).

The geometric path is a bit-faithful reimplementation of R3F ``rrf_times.R``;
the geometric tests below freeze that behaviour as the reference contract, and
``tests/test_regression_golden.py`` remains the authority on numerical
agreement.
"""

from __future__ import annotations

import inspect
import logging
import math
from typing import Callable, Optional

import pytest

from openreltime import _constants as C
from openreltime.rates import rrf_rates
from openreltime.rrf import RrfEngine, RrfRow
from openreltime.times import rrf_rates_times, rrf_times
from openreltime.tree import PhyloNode
from openreltime.treeio import parse_newick

TOL = 1e-12

#: the two averaging conventions of msy044, transcribed from the paper:
#: geometric ``sqrt(b1*b2)`` (eqs 28-42, R3F) vs arithmetic ``1/2*(b1+b2)``
#: (eqs 1-27; fig. 2 legend: "La = b5 + 1/2(b1 + b2) ... arithmetic mean").
FOLD: dict[str, Callable[[float, float], float]] = {
    "geometric": lambda a, b: math.sqrt(a * b),
    "arithmetic": lambda a, b: (a + b) / 2.0,
}

# three-level tree: root -> (node with a cherry + one tip) / cherry.  The left
# child's own composite is sqrt(0.223205*0.20)+0.10 = 0.311284 geometric and
# (0.223205+0.20)/2+0.10 = 0.311603 arithmetic; the two conventions agree only
# at depth 1, so a two-level tree cannot tell them apart.
THREE_LEVEL = "(((A:0.223205,B:0.20):0.10,C:0.15):0.05,(D:0.18,E:0.22):0.12):0.0;"
FOLD_GEO_LEFT = 0.31128416883429766
FOLD_ARI_LEFT = 0.3116025

# three-level tree whose two level-2 nodes each hold two cherries, so no
# tip-grandchild backfill disturbs their local solution and b = r * t can be
# checked at every internal node.
BALANCED3 = (
    "(((A:0.12,B:0.20):0.10,(C:0.05,D:0.30):0.04):0.05,"
    "((E:0.18,F:0.22):0.12,(G:0.09,H:0.31):0.07):0.03);"
)
BALANCED3_LEFT_GEOM = 0.2549193338482967  # sqrt(0.12*0.20) + 0.10
BALANCED3_LEFT_ARI = 0.26  # (0.12+0.20)/2 + 0.10

# two-level tree (four tips, children are cherries): the only depth at which
# the two foldings coincide, so it cannot discriminate between the conventions.
SHALLOW = "((T1:0.4,T2:0.2):0.25,(T3:0.5,T4:0.3):0.35):0.0;"

# six-tip tree that DOES trigger the guard at the default threshold: the
# adjusted rate of the A/B cherry's parent row is 36.77 > 20.
GUARD_TREE = (
    "((((A:0.30,B:0.02):0.30,C:0.02):0.30,D:0.02):0.30,(E:0.05,F:0.06):0.05):0.0;"
)
#: cherry whose age the guard replaces (its unguarded age, its guarded age, and
#: the ages of the three candidate ancestors, all mutually distinct)
GUARD_CHERRY_RAW = 0.0021067162135748638  # t7.adjust without the guard
GUARD_CHERRY_KEPT = 0.010265993023935459  # == parent's t7.adjust, the supplier
GUARD_GRANDPARENT = 0.02375532978011225  # grandparent: must NOT be picked
GUARD_ROOT = 0.20161323649140578  # root: must NOT be picked

# pectinate eight-tip tree: three chained replacements, all supplied by the
# nearest *clean* ancestor (the cherry-line node 10), not by the root.
CHAIN_TREE = (
    "((((((T1:1.11713,T2:0.0006594),T3:0.000801268),T4:0.00132357),"
    "T5:0.0972749),T6:1.09272),(R1:0.805467,R2:1.31768)):0.0;"
)
CHAIN_SUPPLIER_TIME = 1.0302173346241075
CHAIN_ROOT_TIME = 0.3664331954752665


def _engine(newick: str, mean: str, *, threshold: Optional[float] = None) -> RrfEngine:
    tree = parse_newick(newick)
    engine = RrfEngine(
        tree, mean=mean, rate_ratio_threshold=threshold, compute_times=True
    )
    engine.run()
    return engine


def _node(tree: PhyloNode, *tips: str) -> PhyloNode:
    return tree.mrca(list(tips))


# ---------------------------------------------------------------------------
# 07 -- the folding operator must follow the chosen averaging convention
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mean,expected", [("geometric", FOLD_GEO_LEFT),
                                           ("arithmetic", FOLD_ARI_LEFT)])
def test_three_level_composite_follows_mean(mean: str, expected: float) -> None:
    engine = _engine(THREE_LEVEL, mean)
    root_row = engine.row(engine.tree)
    child_row = engine.row(root_row.child1)
    # the parent's composite length is the child's own folded lineage length
    assert root_row.l1 == pytest.approx(expected, abs=TOL)
    assert root_row.l1 == pytest.approx(
        FOLD[mean](child_row.l1, child_row.l2) + child_row.l5, abs=TOL
    )
    # the two conventions disagree at depth 3 by ~4e-4, i.e. the frozen
    # expectations are not degenerate
    assert abs(FOLD_GEO_LEFT - FOLD_ARI_LEFT) > 1e-5


def test_arithmetic_identity_b_equals_r_times_t() -> None:
    """RRF identity ``b = r * t`` for the arithmetic convention.

    Every internal row's composite lineage length must equal rate x age in the
    same local frame: a row folded arithmetically has to be solved
    arithmetically, or the two sides diverge (0.2549193338482967 vs 0.26 at
    the root of :data:`BALANCED3`).
    """
    engine = _engine(BALANCED3, "arithmetic")
    assert engine.row(engine.tree).l1 == pytest.approx(BALANCED3_LEFT_ARI, abs=TOL)
    for row in engine.rows.values():
        la = FOLD["arithmetic"](row.l1, row.l2) + row.l5
        lb = FOLD["arithmetic"](row.l3, row.l4) + row.l6
        assert row.r5 is not None and row.r6 is not None and row.t7 is not None
        assert row.r5 * row.t7 == pytest.approx(la, abs=TOL)
        assert row.r6 * row.t7 == pytest.approx(lb, abs=TOL)


def test_geometric_identity_and_fold_are_unchanged() -> None:
    """Guard rail for review section 7: geometric behaviour must not move."""
    for newick, expected in (
        (THREE_LEVEL, FOLD_GEO_LEFT),
        (BALANCED3, BALANCED3_LEFT_GEOM),
    ):
        engine = _engine(newick, "geometric")
        assert engine.row(engine.tree).l1 == pytest.approx(expected, abs=TOL)
    engine = _engine(BALANCED3, "geometric")
    for row in engine.rows.values():
        la = FOLD["geometric"](row.l1, row.l2) + row.l5
        lb = FOLD["geometric"](row.l3, row.l4) + row.l6
        assert row.r5 is not None and row.r6 is not None and row.t7 is not None
        assert row.r5 * row.t7 == pytest.approx(la, abs=TOL)
        assert row.r6 * row.t7 == pytest.approx(lb, abs=TOL)
    # geometric balances sister rate changes exactly: r5 * r6 = 1 for the rows
    # that still carry their own local rates (tip-grandchild backfills store
    # great-grandchild rates, whose product is not 1)
    for row in engine.rows.values():
        if row.backfilled:
            continue
        assert row.r5 * row.r6 == pytest.approx(1.0, abs=TOL)


def test_two_level_tree_folding_is_convention_independent() -> None:
    """At depth 2 both conventions fold to the tip lengths."""
    geo_engine = _engine(SHALLOW, "geometric")
    ari_engine = _engine(SHALLOW, "arithmetic")
    geo_row = geo_engine.row(geo_engine.tree)
    ari_row = ari_engine.row(ari_engine.tree)
    assert [geo_row.l1, geo_row.l2, geo_row.l3, geo_row.l4] == [
        ari_row.l1,
        ari_row.l2,
        ari_row.l3,
        ari_row.l4,
    ]
    assert geo_row.l1 == pytest.approx(0.4, abs=TOL)
    assert geo_row.l2 == pytest.approx(0.2, abs=TOL)


# ---------------------------------------------------------------------------
# 30 -- one folding helper, used by every call site
# ---------------------------------------------------------------------------

def test_all_folding_sites_route_through_composite(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """``_composite`` is the only folding path; no site folds inline."""
    calls: list[int] = []
    original = RrfEngine._composite

    def spy(self: RrfEngine, row: RrfRow) -> tuple[float, float]:
        calls.append(row.node_id)
        return original(self, row)

    monkeypatch.setattr(RrfEngine, "_composite", spy)
    _engine(BALANCED3, "geometric")
    # every internal child is folded into its parent exactly once: the root's
    # two children plus the two level-2 nodes' two children each
    assert len(calls) == 6
    assert sorted(calls) == [10, 11, 12, 13, 14, 15]


# ---------------------------------------------------------------------------
# 16 -- the registered default threshold must actually be exercised
# ---------------------------------------------------------------------------

def test_default_threshold_triggers_guard_and_names_the_ancestor() -> None:
    tree = parse_newick(GUARD_TREE)
    result = rrf_times(tree)  # no explicit threshold -> the registered default
    assert result.n_rate_guarded == 1

    cherry = _node(tree, "A", "B")
    parent = _node(tree, "A", "B", "C")
    grandparent = _node(tree, "A", "B", "C", "D")
    assert result.times[cherry.node_id] == pytest.approx(
        GUARD_CHERRY_KEPT, abs=TOL
    )
    # the supplier is the nearest clean ancestor = the parent...
    assert result.times[cherry.node_id] == pytest.approx(
        result.times[parent.node_id], abs=TOL
    )
    # ... and not the grandparent, not the root (the failure mode the
    # membership-only assertion could not see)
    assert result.times[grandparent.node_id] == pytest.approx(
        GUARD_GRANDPARENT, abs=TOL
    )
    assert result.times[cherry.node_id] != pytest.approx(
        result.times[grandparent.node_id], abs=1e-6
    )
    assert result.times[cherry.node_id] != pytest.approx(
        result.times[tree.node_id], abs=1e-6
    )
    assert result.times[tree.node_id] == pytest.approx(GUARD_ROOT, abs=TOL)


def test_guard_disabled_keeps_the_original_age() -> None:
    tree = parse_newick(GUARD_TREE)
    result = rrf_times(tree, rate_ratio_threshold=None)
    cherry = _node(tree, "A", "B")
    assert result.n_rate_guarded == 0
    assert result.times[cherry.node_id] == pytest.approx(GUARD_CHERRY_RAW, abs=TOL)


def test_guard_chain_resolves_to_the_nearest_clean_ancestor() -> None:
    tree = parse_newick(CHAIN_TREE)
    result = rrf_times(tree)
    assert result.n_rate_guarded == 3
    replaced = [_node(tree, "T1", "T2", "T3"), _node(tree, "T1", "T2", "T3", "T4")]
    for node in replaced:
        assert result.times[node.node_id] == pytest.approx(
            CHAIN_SUPPLIER_TIME, abs=TOL
        )
        assert result.times[node.node_id] != pytest.approx(
            result.times[tree.node_id], abs=1e-6
        )
    supplier = _node(tree, "T1", "T2", "T3", "T4", "T5", "T6")
    assert result.times[supplier.node_id] == pytest.approx(
        CHAIN_SUPPLIER_TIME, abs=TOL
    )
    assert result.times[tree.node_id] == pytest.approx(CHAIN_ROOT_TIME, abs=TOL)


def test_guard_rewrites_times_and_never_rates() -> None:
    """R3F rrf_times.R:353-373 touches ``t7.adjust`` only, never the rates."""
    guarded = rrf_times(parse_newick(GUARD_TREE))
    unguarded = rrf_times(parse_newick(GUARD_TREE), rate_ratio_threshold=None)
    assert guarded.n_rate_guarded == 1 and unguarded.n_rate_guarded == 0
    cherry = _node(guarded.tree, "A", "B")
    assert guarded.times[cherry.node_id] != pytest.approx(
        unguarded.times[cherry.node_id], abs=1e-9
    )
    assert set(guarded.rates) == set(unguarded.rates)
    for node_id, rate in guarded.rates.items():
        assert rate == pytest.approx(unguarded.rates[node_id], abs=TOL)


# ---------------------------------------------------------------------------
# 18 -- rrf_rates cannot honour a threshold and must say so
# ---------------------------------------------------------------------------

def test_rrf_rates_ignores_threshold_and_warns(caplog) -> None:
    tree = parse_newick(GUARD_TREE)
    with caplog.at_level(logging.WARNING, logger="openreltime"):
        guarded = rrf_rates(tree, rate_ratio_threshold=C.RATE_RATIO_THRESHOLD)
    plain = rrf_rates(tree)
    # the reference rrf_rates.R has no guard: the rates are identical even
    # though this tree holds a node rate of 142 > 20
    assert guarded.rates == plain.rates
    assert any(n > C.RATE_RATIO_THRESHOLD for n in plain.rates.values())
    assert any("ignored by rrf_rates" in r.message for r in caplog.records)
    assert guarded.warnings and "ignored by rrf_rates" in guarded.warnings[0]
    assert not plain.warnings


# ---------------------------------------------------------------------------
# 33 -- one source of default values
# ---------------------------------------------------------------------------

def test_guard_defaults_reference_the_constants_module() -> None:
    """The 20.0/None guard defaults are read from ``_constants``."""
    sources = {
        "times": inspect.getsource(inspect.getmodule(rrf_times)),
        "rates": inspect.getsource(inspect.getmodule(rrf_rates)),
        "rrf": inspect.getsource(inspect.getmodule(RrfEngine)),
    }
    times_default = "rate_ratio_threshold: Optional[float] = C.RATE_RATIO_THRESHOLD,"
    rates_default = (
        "rate_ratio_threshold: Optional[float] = C.RATE_RATIO_GUARD_DISABLED,"
    )
    for name, needle in (
        ("times", times_default),
        ("rates", rates_default),
        ("rrf", rates_default),
    ):
        assert needle in sources[name], f"{name}.py does not use the constant"
        assert "= 20.0" not in sources[name], f"{name}.py still hard-codes 20.0"
    for func in (rrf_times, rrf_rates_times):
        default = inspect.signature(func).parameters["rate_ratio_threshold"].default
        assert default == C.RATE_RATIO_THRESHOLD and isinstance(default, float)
    for func in (rrf_rates, RrfEngine.__init__):
        params = inspect.signature(func).parameters
        assert params["rate_ratio_threshold"].default is C.RATE_RATIO_GUARD_DISABLED
    assert C.RATE_RATIO_GUARD_DISABLED is None
