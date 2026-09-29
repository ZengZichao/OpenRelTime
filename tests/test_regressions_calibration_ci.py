"""Regression tests for the calibration and confidence-interval layer.

Each test pins one contract the calibration and CI code must satisfy:

* a calibration row carrying neither bound is rejected when the table is
  parsed, naming the offending line, so no infinite time factor can appear;
* ``RV(R)`` in msz236 eq (14) takes its average and its sum on the same
  per-lineage scale, so the rate-heterogeneity term is positive in the
  few-sites regime rather than vanishing identically;
* a node pinned at a hard bound carries the uncertainty its own calibration
  asserts, so its interval is never of exactly zero width;
* the differentiation pipeline follows the ``mean`` convention and the effective
  bounds that produced the point estimate, so derivatives and point estimate
  always come from one estimator;
* the per-lineage rate reported by ``calibrate`` and the RRF rate are distinct
  quantities with distinct names, and ``Vobs(R)`` uses the RRF one;
* the sampling-variance clamp warns when it actually binds;
* the ``times`` branch indexes trial-pipeline output only after checking that it
  carries the analysed tree's own node ids.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pytest

from openreltime.calibrate import (
    Calibration,
    _CalibratedTree,
    calibrate,
    parse_calibrations,
)
from openreltime.ci import _make_pipeline, confidence_interval
from openreltime.times import rrf_times
from openreltime.tree import PhyloNode
from openreltime.treeio import parse_newick

# A 7-taxon ingroup tree with strong rate heterogeneity.  After id assignment
# (tips 1-7, internal nodes cladewise) node 8 is the root, node 9 the
# (A,B,C,D) clade, node 10 the (A,B) clade, node 12 the (E,F,G) clade.
SEVEN_TAXA = (
    "(((A:0.05,B:0.25):0.15,(C:0.30,D:0.04):0.10):0.20,"
    "((E:0.06,F:0.28):0.12,G:0.35):0.08):0.0;"
)
# a near-clocklike tree: its sampling variance swamps the observed rate variance
CLOCKLIKE = (
    "(((A:0.10,B:0.12):0.20,(C:0.15,D:0.18):0.22):0.30,"
    "((E:0.20,F:0.22):0.15,G:0.28):0.25):0.0;"
)
# node 10's time is replaced by its ancestor's by the rate-ratio guard, so its
# calibrated duration is exactly zero
GUARD_TREE = (
    "(((A:1e-6,B:0.9):0.001,(C:0.4,D:0.5):0.3):0.1,"
    "(E:0.2,(F:0.1,G:0.3):0.2):0.4):0.0;"
)


def _tree(newick: str) -> PhyloNode:
    tree = parse_newick(newick)
    tree.assign_ids()
    return tree


def _blens(tree: PhyloNode) -> np.ndarray:
    return np.array(
        [float(n.blen or 0.0) for n in tree.walk() if not n.is_root()], dtype=float
    )


def _pinned_case():
    """7-taxon tree whose two calibrated nodes are both pinned at a hard bound.

    The bound intervals of the ancestor (node 9) and the descendant (node 10)
    are mutually incompatible with the tree's relative-time ratio, so the
    violation loop drives node 9 onto its maximum and node 10 onto its minimum.
    """
    tree = _tree(SEVEN_TAXA)
    times = rrf_times(tree)
    cals = [
        Calibration(node_id=9, min_bound=0.4, max_bound=0.5),
        Calibration(node_id=10, min_bound=0.3, max_bound=0.35),
    ]
    cal = calibrate(times, cals, method="bounds")
    return tree, times, cal


def _row(ci, node_id: int):
    return ci.table[ci.table["node_id"] == node_id].iloc[0]


# ---------------------------------------------------------------------------
#  - unbounded calibration rows
# ---------------------------------------------------------------------------

def test_calibration_without_any_bound_is_rejected() -> None:
    """A row with neither min nor max must fail at construction, not at solve."""
    with pytest.raises(ValueError, match="at least one finite bound"):
        Calibration(node_id=9)
    with pytest.raises(ValueError, match="at least one finite bound"):
        Calibration(taxon_set=frozenset({"A", "B"}))
    # non-finite and inverted intervals are rejected too
    with pytest.raises(ValueError, match="finite"):
        Calibration(node_id=9, min_bound=float("inf"))
    with pytest.raises(ValueError, match="max_bound must be"):
        Calibration(node_id=9, min_bound=0.0, max_bound=0.0)
    with pytest.raises(ValueError, match="interval is empty"):
        Calibration(node_id=9, min_bound=0.9, max_bound=0.4)


def test_parse_calibrations_names_the_offending_line(tmp_path: Path) -> None:
    path = tmp_path / "calibrations.tsv"
    path.write_text(
        "node_id\ttaxon_set\tmin_bound\tmax_bound\tdensity\tdensity_params\n"
        "\tA|B\t0.4\t0.5\t.\t.\n"
        "\tC|D\t.\t.\t.\t.\n"
    )
    with pytest.raises(ValueError, match=r"line 3.*at least one finite bound"):
        parse_calibrations(path)


def test_bounds_method_rejects_density_only_rows() -> None:
    """A density-only row with no finite bound must not yield f = inf."""
    tree = _tree(SEVEN_TAXA)
    times = rrf_times(tree)
    cals = [
        Calibration(
            node_id=9, density="uniform", density_params={"min": 0.4, "max": 0.5}
        )
    ]
    with pytest.raises(ValueError, match=r"method='bounds' cannot use density-only"):
        calibrate(times, cals, method="bounds")
    # the same row is legitimate under the effective method (msz236)
    cal = calibrate(times, cals, method="effective", n_effective=40, seed=7)
    assert np.isfinite(cal.time_factor)
    assert all(np.isfinite(v) for v in cal.times.values())


def test_feasible_factor_raises_instead_of_returning_inf() -> None:
    """Depth defence: the solver itself refuses to produce a non-finite factor."""
    tree = _tree(SEVEN_TAXA)
    times = rrf_times(tree)
    node = next(n for n in tree.walk() if n.node_id == 9)
    solver = _CalibratedTree(tree, times, [(node, None, None)])
    with pytest.raises(ValueError, match=r"no finite calibration bound for node\(s\) \[9\]"):
        solver.solve()


# ---------------------------------------------------------------------------
# the rate-heterogeneity component must not be scaled away
# ---------------------------------------------------------------------------

def test_rate_heterogeneity_variance_is_positive_at_low_site_counts() -> None:
    """``RV(R)`` stays positive at a low site count (7 taxa, L = 30)."""
    tree, _times, cal = _pinned_case()
    ci = confidence_interval(cal, seq_length=30)
    assert ci.params["n_lineages"] == 12
    # Vobs(R) of eq (9) is an average over lineages, SV(R) of eq (11) a sum: the
    # subtraction happens per lineage, so RV(R) survives even at L=30.
    assert ci.params["Vobs_R"] > 0.0
    assert 0.0 < ci.params["RV_R"] < ci.params["Vobs_R"]
    assert ci.params["SV_R"] == pytest.approx(
        ci.params["SV_per_lineage_R"] * ci.params["n_lineages"]
    )
    assert ci.params["SV_R"] > ci.params["Vobs_R"]  # eq (11) sums, eq (9) averages
    # vR(b) of eq (14) therefore contributes: the branch variances are strictly
    # larger than the pure sampling variances b/L
    blens = _blens(tree)
    v_s = blens / 30.0
    v_r_total = ci.params["RV_R"]  # sum_j vR(b_j) = RV(R) by eq (14)
    assert v_r_total > 0.0
    assert float(np.sum(v_s) + v_r_total) > float(np.sum(v_s))


def test_rate_heterogeneity_saturates_at_Vobs_as_sites_grow() -> None:
    """eq (13): SV(R)/N ~ 1/L, so RV(R) rises monotonically towards Vobs(R)."""
    _tree_, _times, cal = _pinned_case()
    rvs = []
    for length in (30, 100, 300, 1000, 10000):
        ci = confidence_interval(cal, seq_length=length)
        rvs.append(ci.params["RV_R"])
        assert ci.params["Vobs_R"] == pytest.approx(rvs[-1] + ci.params["SV_per_lineage_R"])
    assert all(b > a for a, b in zip(rvs, rvs[1:])), rvs
    assert rvs[0] > 0.0
    assert rvs[-1] < confidence_interval(cal).params["Vobs_R"]
    # in the sampling-dominated regime (few sites) RV(R) is clamped away and the
    # interval is set by vS(b) = b/L, so it does widen as sites are removed
    ci_few = confidence_interval(cal, seq_length=5)
    ci_noiseless = confidence_interval(cal)
    assert ci_few.params["RV_R"] == 0.0
    assert float(ci_few.table["width"].mean()) > float(ci_noiseless.table["width"].mean())


# ---------------------------------------------------------------------------
#  - hard-bound-pinned nodes
# ---------------------------------------------------------------------------

def test_pinned_nodes_get_the_calibration_uncertainty() -> None:
    tree, _times, cal = _pinned_case()
    assert cal.times[9] == pytest.approx(0.5)  # pinned at its maximum bound
    assert cal.times[10] == pytest.approx(0.3)  # pinned at its minimum bound
    ci = confidence_interval(cal, seq_length=30)
    for node_id, lo, hi in ((9, 0.4, 0.5), (10, 0.3, 0.35)):
        row = _row(ci, node_id)
        assert row["width"] > 0.0, "a calibrated node must not report zero uncertainty"
        assert row["lower"] == pytest.approx(lo)
        assert row["upper"] == pytest.approx(hi)
        assert row["se"] > 0.0
        assert not bool(row["se_reliable"])
        assert "pinned" in row["notes"]
    # nodes qualified by eq (7) are untouched and flagged as reliable
    root_row = _row(ci, [n.node_id for n in tree.walk() if n.is_root()][0])
    assert bool(root_row["se_reliable"])
    assert root_row["notes"] == ""
    # every interval contains its point estimate and is non-negative
    assert (ci.table["lower"] <= ci.table["time"] + 1e-12).all()
    assert (ci.table["upper"] >= ci.table["time"] - 1e-12).all()
    assert (ci.table["lower"] >= 0.0).all()
    assert (ci.table["width"] >= 0.0).all()


def test_midpoint_fixed_calibration_is_not_zero_width() -> None:
    """The common case: one min-max row fixes its node at the interval midpoint.

    The age is then independent of every branch length, so eq (7) gives exactly
    ``se = 0``; the node must inherit the calibration interval instead of being
    printed with zero uncertainty.
    """
    tree = _tree(CLOCKLIKE)
    times = rrf_times(tree, normalize=True)
    target = tree.mrca(["A", "B"])
    cal = calibrate(
        times,
        [Calibration(node_id=target.node_id, min_bound=1.0, max_bound=2.0)],
        method="bounds",
    )
    assert cal.times[target.node_id] == pytest.approx(1.5)
    ci = confidence_interval(cal, seq_length=500)
    row = _row(ci, target.node_id)
    assert row["width"] == pytest.approx(1.0)
    assert (row["lower"], row["upper"]) == (pytest.approx(1.0), pytest.approx(2.0))
    assert not bool(row["se_reliable"])
    assert "calibration interval" in row["notes"]


def test_one_sided_bound_keeps_truncation_but_not_zero_width() -> None:
    """A min-only calibration pins its node; the open side must stay positive."""
    tree = _tree(CLOCKLIKE)
    times = rrf_times(tree)
    cal = calibrate(
        times, [Calibration(node_id=9, min_bound=0.5, max_bound=None)], method="bounds"
    )
    assert cal.times[9] == pytest.approx(0.5)  # feasible factor sits on the bound
    ci = confidence_interval(cal, seq_length=1000)
    row = _row(ci, 9)
    assert row["lower"] == pytest.approx(0.5)  # truncation at the hard bound kept
    assert row["width"] > 0.0
    assert not bool(row["se_reliable"])
    assert "one-sided hard bound" in row["notes"]


def test_unsatisfied_constraint_never_inverts_the_interval() -> None:
    """A bound the solver cannot reach must not print a negative width."""
    tree = _tree(CLOCKLIKE)
    times = rrf_times(tree)
    cal = calibrate(
        times,
        [
            Calibration(node_id=10, min_bound=0.5, max_bound=0.9),
            Calibration(node_id=12, min_bound=0.1, max_bound=0.3),
        ],
        method="bounds",
    )
    assert cal.warnings, "the solver must report the constraint it could not meet"
    ci = confidence_interval(cal, seq_length=30)
    assert (ci.table["width"] >= 0.0).all()
    assert (ci.table["lower"] <= ci.table["time"] + 1e-12).all()
    assert (ci.table["upper"] >= ci.table["time"] - 1e-12).all()


# ---------------------------------------------------------------------------
#  - the derivatives must come from the point estimate's own estimator
# ---------------------------------------------------------------------------

def test_arithmetic_point_estimate_gets_arithmetic_derivatives() -> None:
    tree = _tree(SEVEN_TAXA)
    cals = [Calibration(node_id=9, min_bound=0.4, max_bound=0.9)]
    geo = calibrate(rrf_times(tree), cals, method="bounds")
    ari = calibrate(
        rrf_times(tree, mean="arithmetic"), cals, method="bounds"
    )
    assert ari.params["mean"] == "arithmetic"
    blens = _blens(tree)
    internal = [n.node_id for n in tree.walk() if not n.is_tip()]
    for cal in (geo, ari):
        pipeline = _make_pipeline(cal, "times")
        # the differentiation pipeline reproduces the reported ages exactly, i.e.
        # it re-runs the same estimator rather than the geometric default
        assert np.allclose(
            pipeline(blens), [cal.times[nid] for nid in internal], rtol=1e-9
        ), cal.params["mean"]
    ci_geo = confidence_interval(geo, seq_length=200)
    ci_ari = confidence_interval(ari, seq_length=200)
    assert ci_ari.params["pipeline_mean"] == "arithmetic"
    assert ci_geo.params["pipeline_mean"] == "geometric"
    assert not np.allclose(
        ci_geo.table["se"].to_numpy(), ci_ari.table["se"].to_numpy(), rtol=1e-6
    ), "arithmetic and geometric estimates must not share one derivative path"


def test_effective_result_is_differentiated_along_its_effective_bounds() -> None:
    tree = _tree(SEVEN_TAXA)
    times = rrf_times(tree)
    cals = [
        Calibration(
            node_id=9, density="uniform", density_params={"min": 0.4, "max": 0.9}
        )
    ]
    cal = calibrate(times, cals, method="effective", n_effective=200, seed=11)
    assert cal.effective_bounds is not None
    blens = _blens(tree)
    internal = [n.node_id for n in tree.walk() if not n.is_tip()]
    pipeline = _make_pipeline(cal, "times")
    assert np.allclose(
        pipeline(blens), [cal.times[nid] for nid in internal], rtol=1e-9
    ), "the CI pipeline must re-run the effective bounds, not the raw rows"
    ci = confidence_interval(cal, seq_length=500)
    assert ci.params["pipeline_method"] == "effective"
    assert np.all(np.isfinite(ci.table["time"]))
    assert (ci.table["width"] > 0.0).all()


def test_unknown_estimator_is_reported() -> None:
    _tree_, _times, cal = _pinned_case()
    cal.params["method"] = "bootstrap"
    with pytest.raises(ValueError, match="method='bounds' or method='effective'"):
        confidence_interval(cal, seq_length=30)
    cal.params["method"] = "bounds"
    cal.params["mean"] = "harmonic"
    with pytest.raises(ValueError, match="averaging convention"):
        confidence_interval(cal, seq_length=30)


# ---------------------------------------------------------------------------
#  - two different rate vectors, one of them used for the variance
# ---------------------------------------------------------------------------

def test_rrf_and_implied_rates_are_named_apart() -> None:
    tree, times, cal = _pinned_case()
    frame = cal.to_pandas()
    assert {"RRFRate", "ImpliedRate"} <= set(frame.columns)
    assert "Rate" not in frame.columns
    # rates is the RRF quantity of eqs (1)-(4), identical to rrf_times' Rate column
    assert cal.rates == pytest.approx({k: v for k, v in times.rates.items()})
    # the implied absolute rate is a different number for internal lineages
    internal_edges = [n for n in tree.walk() if not n.is_root() and not n.is_tip()]
    assert any(
        abs(cal.rates[n.node_id] - cal.implied_rates[n.node_id]) > 1e-6
        for n in internal_edges
    )


def test_variance_uses_the_rrf_rates_and_drops_lineages_without_duration() -> None:
    tree = _tree(GUARD_TREE)
    times = rrf_times(tree)
    assert times.n_rate_guarded > 0, "the fixture must produce a zero-duration lineage"
    zero_nodes = [
        n.node_id
        for n in tree.walk()
        if not n.is_root() and times.times[n.parent.node_id] <= times.times[n.node_id]
    ]
    assert zero_nodes
    cal = calibrate(
        times,
        [Calibration(node_id=tree.node_id, min_bound=1.0, max_bound=2.0)],
        method="bounds",
    )
    for nid in zero_nodes:
        assert np.isnan(cal.implied_rates[nid]), "no 0.0 placeholder for duration <= 0"
    assert cal.params["n_lineages_dropped"] == len(zero_nodes)
    assert any("non-positive" in w for w in cal.warnings)

    ci = confidence_interval(cal, seq_length=100)
    edge_nodes = [n for n in tree.walk() if not n.is_root()]
    rrf = np.array([cal.rates[n.node_id] for n in edge_nodes], dtype=float)
    assert ci.params["rate_definition"].startswith("RRF")
    assert ci.params["Vobs_R"] == pytest.approx(
        float(np.mean((rrf - rrf.mean()) ** 2)), rel=1e-9
    )
    assert ci.params["n_lineages"] == ci.params["n_branches"] - ci.params[
        "n_lineages_dropped"
    ]
    assert np.isfinite(ci.params["Vobs_R"])


def test_non_finite_rates_are_dropped_from_the_average() -> None:
    _tree_, _times, cal = _pinned_case()
    victim = sorted(cal.rates)[3]
    cal.rates[victim] = float("nan")
    ci = confidence_interval(cal, seq_length=300)
    assert ci.params["n_lineages_dropped"] == 1
    assert ci.params["n_lineages"] == ci.params["n_branches"] - 1
    assert np.isfinite(ci.params["Vobs_R"]) and np.isfinite(ci.params["RV_R"])


# ---------------------------------------------------------------------------
# the clamp warning must be reachable at all
# ---------------------------------------------------------------------------

def test_clamp_warning_fires_when_sampling_variance_exceeds_vobs(
    caplog: pytest.LogCaptureFixture,
) -> None:
    tree = _tree(CLOCKLIKE)
    times = rrf_times(tree)
    cal = calibrate(
        times, [Calibration(node_id=9, min_bound=0.4, max_bound=0.9)], method="bounds"
    )
    with caplog.at_level(logging.WARNING, logger="openreltime"):
        ci = confidence_interval(cal, seq_length=10)
    assert ci.params["RV_R"] == 0.0
    assert ci.params["SV_per_lineage_R"] > ci.params["Vobs_R"]
    messages = [r.getMessage() for r in caplog.records]
    assert any(
        "sampling variance exceeds observed rate variance" in m for m in messages
    ), messages
    assert any("clamped to 0" in m for m in messages)


# ---------------------------------------------------------------------------
#  - explicit node-id alignment in the pipeline
# ---------------------------------------------------------------------------

def test_divergent_trial_ids_are_detected(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _tree_, _times, cal = _pinned_case()
    tree = _tree(SEVEN_TAXA)
    blens = _blens(tree)
    pipeline = _make_pipeline(cal, "times")
    assert len(pipeline(blens)) == len([n for n in tree.walk() if not n.is_tip()])

    def scrambled(self: PhyloNode) -> None:
        for i, node in enumerate(reversed(list(self.walk()))):
            node.node_id = i

    monkeypatch.setattr(PhyloNode, "assign_ids", scrambled)
    with pytest.raises(ValueError, match="node-id space"):
        pipeline(blens)


# ---------------------------------------------------------------------------
# provenance hand-off: unrecognised keys must be absorbed and echoed back
# ---------------------------------------------------------------------------

def test_provenance_keys_flow_from_calibrate_into_the_ci() -> None:
    tree = _tree(SEVEN_TAXA)
    times = rrf_times(tree)
    provenance = {
        "tree_file": "/tmp/seven.nwk",
        "calibrations_file": "/tmp/seven_calibrations.tsv",
        "outgroup": ["X", "Y"],
        "input_fmt": "nexus",
        "resolve_polytomy": "random",
        "outgroup_check": "warn",
    }
    cal = calibrate(
        times,
        [Calibration(node_id=9, min_bound=0.4, max_bound=0.9)],
        method="bounds",
        **provenance,
    )
    for key, value in provenance.items():
        assert cal.params[key] == value
    ci = confidence_interval(cal, seq_length=500)
    assert ci.params["outgroup_check"] == "warn"
    ci2 = confidence_interval(cal, seq_length=500, outgroup_check="error")
    assert ci2.params["outgroup_check"] == "error"
    assert ci2.params["level"] == 0.95
