"""Calibration module behaviour, the CI module, and edge cases (validation gate G7)."""

from __future__ import annotations

import logging
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from openreltime.calibrate import (
    Calibration,
    calibrate,
    density_quantiles,
    parse_calibrations,
)
from openreltime.ci import confidence_interval
from openreltime.table import tree2table
from openreltime.times import rrf_times
from openreltime.treeio import parse_newick, read_tree, write_nexus

ROOT = Path(__file__).resolve().parents[1]
TREE = ROOT / "data" / "examples" / "example.nwk"
OUTGROUP = ["Ornithorhynchus_anatinus", "Zaglossus_bruijni", "Tachyglossus_aculeatus"]

SIMPLE = "((A:0.10,B:0.12):0.20,(C:0.15,(D:0.18,E:0.22):0.25):0.30):0.0;"


def _simple_tree():
    return parse_newick(SIMPLE)


# ---------------------------------------------------------------------------
# G7: calibration
# ---------------------------------------------------------------------------

@pytest.mark.gate("G7")
def test_bounds_calibration_hits_bounds() -> None:
    """A single min-max calibration must be satisfied exactly (midpoint f)."""
    tree = _simple_tree()
    times = rrf_times(tree, normalize=True)
    target = tree.mrca(["A", "B"])
    cal = calibrate(
        times,
        [Calibration(node_id=target.node_id, min_bound=1.0, max_bound=2.0)],
        method="bounds",
    )
    age = cal.times[target.node_id]
    assert 1.0 <= age <= 2.0
    # midpoint semantics: f = mid of [1/t, 2/t] -> age = 1.5
    assert age == pytest.approx(1.5, rel=1e-9)
    # all absolute times scale consistently
    assert cal.time_factor > 0


@pytest.mark.gate("G7")
def test_two_calibrations_jointly_satisfied() -> None:
    tree = _simple_tree()
    times = rrf_times(tree, normalize=True)
    n1 = tree.mrca(["A", "B"])
    n2 = tree.mrca(["D", "E"])
    cal = calibrate(
        times,
        [
            Calibration(node_id=n1.node_id, min_bound=0.8, max_bound=1.2),
            Calibration(node_id=n2.node_id, min_bound=0.3, max_bound=0.6),
        ],
        method="bounds",
    )
    for node, (mn, mx) in ((n1, (0.8, 1.2)), (n2, (0.3, 0.6))):
        assert mn <= cal.times[node.node_id] <= mx


@pytest.mark.gate("G7")
def test_conflicting_calibration_raises() -> None:
    """Impossible constraints (child older than parent bounds) must error."""
    tree = _simple_tree()
    times = rrf_times(tree, normalize=True)
    n_child = tree.mrca(["A", "B"])
    n_root = tree  # root of the whole tree
    # root forced young, child forced older than root can be
    with pytest.raises(ValueError):
        calibrate(
            times,
            [
                Calibration(node_id=n_root.node_id, min_bound=None, max_bound=0.5),
                Calibration(node_id=n_child.node_id, min_bound=0.5, max_bound=None),
            ],
            method="bounds",
        )


@pytest.mark.gate("G7")
def test_calibration_by_taxon_set_and_effective(tmp_path: Path) -> None:
    cal_file = tmp_path / "calibrations.tsv"
    cal_file.write_text(
        "node_id\ttaxon_set\tmin_bound\tmax_bound\tdensity\tdensity_params\n"
        "\tA|B\t1.0\t2.0\t\t\n"
        "\tD|E\t.\t.\texponential\toffset=0.3;mean=0.5\n"
    )
    cals = parse_calibrations(cal_file)
    assert len(cals) == 2
    assert cals[0].taxon_set == frozenset(["A", "B"])
    assert cals[1].density_params == {"offset": 0.3, "mean": 0.5}

    tree = _simple_tree()
    times = rrf_times(tree, normalize=True)
    cal = calibrate(times, cals, method="effective", n_effective=500, seed=1)
    assert cal.effective_bounds is not None
    assert len(cal.effective_bounds) == 2
    assert (cal.effective_bounds["effective_min"] <= cal.effective_bounds["effective_max"]).all()
    assert cal.empirical_ci is not None
    n1 = tree.mrca(["A", "B"])
    row = cal.effective_bounds[cal.effective_bounds["node_id"] == n1.node_id]
    assert 1.0 - 1e-9 <= float(row["effective_min"].iloc[0])
    assert float(row["effective_max"].iloc[0]) <= 2.0 + 1e-9

    # the density row must be contained in *its own* parent density, not merely
    # in the hard-bound control row above: the effective bounds are the 2.5/97.5
    # percentiles of the replicate node ages, and those replicate pairs are drawn
    # from the exponential, so the pair can never leave the density's own
    # 2.5/97.5 quantiles (this is the assertion `density_quantiles` existed for)
    n2 = tree.mrca(["D", "E"])
    lo, hi = density_quantiles("exponential", {"offset": 0.3, "mean": 0.5})
    drow = cal.effective_bounds[cal.effective_bounds["node_id"] == n2.node_id]
    eff_min = float(drow["effective_min"].iloc[0])
    eff_max = float(drow["effective_max"].iloc[0])
    assert lo - 1e-9 <= eff_min <= eff_max <= hi + 1e-9, (
        f"effective bounds ({eff_min}, {eff_max}) left the exponential's own "
        f"2.5/97.5 quantiles ({lo}, {hi})"
    )


# ---------------------------------------------------------------------------
# CI module
# ---------------------------------------------------------------------------

@pytest.mark.gate("G7")
def test_ci_identity_for_zero_variance() -> None:
    """vS = 0 and no rate variance collapses the CI onto the point estimate."""
    tree = _simple_tree()
    times = rrf_times(tree, normalize=True)
    target = tree.mrca(["A", "B"])
    # constant-rate tree: all branch lengths equal -> rates equal
    cal = calibrate(
        times,
        [Calibration(node_id=target.node_id, min_bound=1.0, max_bound=2.0)],
        method="bounds",
    )
    ci = confidence_interval(cal, level=0.95)
    row = ci.table[ci.table["node_id"] == target.node_id].iloc[0]
    assert row["lower"] <= row["time"] <= row["upper"]
    assert (ci.table["width"] >= 0).all()


@pytest.mark.gate("G7")
def test_variance_components_move_in_opposite_directions_with_L() -> None:
    """msz236 splits the observed rate variance between sampling and rate shift.

    The interval is not monotone in the sequence length: with ``Vobs(R)`` taken
    from the fitted rates, decreasing the sequence length moves variance *from*
    RV(R) *into* vS(b) rather than adding to it, so the total can fall.  What
    must hold is the decomposition itself.
    """
    tree = _simple_tree()
    times = rrf_times(tree, normalize=True)
    target = tree.mrca(["A", "B"])
    cal = calibrate(
        times,
        [Calibration(node_id=target.node_id, min_bound=1.0, max_bound=2.0)],
        method="bounds",
    )
    long_run = confidence_interval(cal, seq_length=100000)
    short_run = confidence_interval(cal, seq_length=50)
    zero = confidence_interval(cal)

    # no sampling information at all: everything is attributed to rate shift
    assert zero.params["RV_R"] == pytest.approx(zero.params["Vobs_R"], rel=1e-9)
    # fewer sites -> more sampling variance and less recoverable rate variance
    assert short_run.params["RV_R"] < long_run.params["RV_R"]
    assert long_run.params["RV_R"] <= zero.params["RV_R"] + 1e-12
    # and the sampling term itself is monotone in 1/L
    assert (short_run.params["SV_R"] or 0) > (long_run.params["SV_R"] or 0)
    for ci in (zero, long_run, short_run):
        assert (ci.table["width"] > 0).all()
        assert (ci.table["se"] >= 0).all()


@pytest.mark.gate("G7")
def test_ci_truncation_at_hard_bounds() -> None:
    tree = _simple_tree()
    times = rrf_times(tree, normalize=True)
    target = tree.mrca(["A", "B"])
    cal = calibrate(
        times,
        [Calibration(node_id=target.node_id, min_bound=1.0, max_bound=2.0)],
        method="bounds",
    )
    ci = confidence_interval(cal, seq_length=100, level=0.999)
    row = ci.table[ci.table["node_id"] == target.node_id].iloc[0]
    assert row["lower"] >= 1.0 - 1e-9
    assert row["upper"] <= 2.0 + 1e-9


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

@pytest.mark.gate("G4")
def test_negative_branch_length_rejected(tmp_path: Path) -> None:
    bad = tmp_path / "neg.nwk"
    bad.write_text("(A:1,B:-2);")
    with pytest.raises(ValueError, match="negative branch length"):
        read_tree(bad)


@pytest.mark.gate("G4")
def test_duplicate_tip_names_rejected(tmp_path: Path) -> None:
    bad = tmp_path / "dup.nwk"
    bad.write_text("((A:1,A:2):0.5,(C:1,D:2):0.4):0.0;")
    with pytest.raises(ValueError, match="duplicate tip names"):
        read_tree(bad)


@pytest.mark.gate("G4")
def test_polytomy_rejected_and_random_resolution(tmp_path: Path) -> None:
    poly = tmp_path / "poly.nwk"
    poly.write_text("((A:1,B:2,C:3):4,(D:5,E:6):7);")
    with pytest.raises(ValueError, match="polytom"):
        read_tree(poly)
    tree = read_tree(poly, resolve_polytomy="random", seed=1)
    assert tree.is_binary()


@pytest.mark.gate("G4")
def test_two_tip_ingroup_rejected() -> None:
    with pytest.raises(ValueError, match="at least 3"):
        tree = parse_newick("((A:1,B:2),(C:3,D:4)):0.0;")
        from openreltime.treeio import root_outgroup

        root_outgroup(tree, ["C", "D"])


@pytest.mark.gate("G1")
def test_chain_tree_all_three_clades() -> None:
    """A fully pectinate (caterpillar) tree exercises every 3-clade path."""
    n = 8
    # build ((...(A,B1),B2),B3...)
    text = "A:1.0"
    for i in range(1, n - 1):
        text = f"({text},B{i}:{1.0 + i * 0.1})"
    tree = parse_newick(text + ";")
    times = rrf_times(tree)
    assert all(np.isfinite(v) for v in times.times.values())


@pytest.mark.gate("G1")
def test_extreme_rates_trigger_guard_and_no_inf() -> None:
    newick = "((A:0.0000001,B:2.0):0.1,(C:0.5,D:0.6):0.9):0.05;"
    tree = parse_newick(newick)
    res = rrf_times(tree)
    vals = np.array(list(res.times.values()), dtype=float)
    assert np.isfinite(vals).all()


@pytest.mark.gate("G2c")
def test_nexus_roundtrip_with_rate_annotations(tmp_path: Path) -> None:
    """G2c: the annotated NEXUS tree round-trips through our own reader."""
    tree = _simple_tree()
    res = rrf_times(tree)
    nexus = write_nexus(tree, rate_map=res.rates, time_map=res.times)
    path = tmp_path / "timetree.nexus"
    path.write_text(nexus)
    reread = read_tree(path, fmt="nexus")
    assert set(reread.tip_labels()) == set(tree.tip_labels())
    # tips keep their [&rate=...] annotations in the emitted text
    text = path.read_text()
    assert "[&rate=" in text


@pytest.mark.gate("G1")
def test_arithmetic_mean_runs_and_differs() -> None:
    tree = _simple_tree()
    geo = rrf_times(tree)
    ari = rrf_times(tree, mean="arithmetic")
    assert any(
        abs(geo.times[k] - ari.times[k]) > 1e-9 for k in geo.times
    )


@pytest.mark.gate("G1")
def test_normalize_option() -> None:
    tree = _simple_tree()
    res = rrf_times(tree, normalize=True)
    assert max(res.times.values()) == pytest.approx(1.0)


def test_tree2table_time_column() -> None:
    tree = _simple_tree()
    frame = tree2table(tree, time=True)
    assert "Time" in frame.columns
    root_row = frame[frame["NodeLabel"] == "-"].iloc[0]
    assert float(root_row["Time"]) > 0


# ---------------------------------------------------------------------------
# G7: the converged solution must still be nested in the tree
# ---------------------------------------------------------------------------

#: The three rows of ``data/examples/example_calibrations.tsv`` are placeholder
#: values that no rate rescaling of the benchmark can satisfy jointly: the
#: ``Elephas``/``Loxodonta`` density implies a far older xenarthran node than
#: that row's own 70-105 bound allows.  Under ``method="effective"`` the repair
#: loop used to converge on such a solution anyway and report it as warnings.
_FORMAT_DEMO = ROOT / "data" / "examples" / "example_calibrations.tsv"


@pytest.mark.gate("G7")
def test_effective_method_rejects_the_non_nested_format_demo() -> None:
    """The documented pathology of the format-illustration file now raises.

    Reproduced before the fix (``n_effective=10000``, ``seed=20260906``):
    ``time_factor`` 4432.8, root 1162.86 Ma, node 307 at -3998.56 Ma, four
    parent-child inversions and four ``NaN`` implied rates, with nothing but
    four log warnings.  ``calibrate()`` must refuse that instead of emitting it.
    """
    tree = read_tree(TREE, outgroup=OUTGROUP)
    times = rrf_times(tree)
    cals = parse_calibrations(_FORMAT_DEMO)
    with pytest.raises(ValueError, match="not nested in the tree") as exc:
        calibrate(times, cals, method="effective", n_effective=200, seed=20260906)
    # the offending nodes are named, not just the fact that something failed
    assert "307" in str(exc.value)


@pytest.mark.gate("G7")
def test_feasible_example_calibrations_still_solve_cleanly() -> None:
    """Rejecting the pathological file must not break the runnable example."""
    tree = read_tree(TREE, outgroup=OUTGROUP)
    times = rrf_times(tree)
    cal = calibrate(
        times, parse_calibrations(ROOT / "data" / "examples" / "example_calibrations_feasible.tsv")
    )
    assert cal.time_factor == pytest.approx(371.59973111298007, rel=1e-12)
    assert min(cal.times.values()) == 0.0
    assert cal.times[tree.node_id] == pytest.approx(97.48139749375025, rel=1e-12)
    assert cal.warnings == []
    assert not any(np.isnan(v) for v in cal.implied_rates.values())
    for node in tree.walk():
        if node.parent is not None:
            assert cal.times[node.node_id] <= cal.times[node.parent.node_id]


@pytest.mark.gate("G7")
def test_time_order_check_names_negative_and_inverted_nodes() -> None:
    """The guard is a check over *every* emitted age, not over one violator."""
    from openreltime.calibrate import _CalibratedTree

    tree = _simple_tree()
    times = rrf_times(tree, normalize=True)
    solver = _CalibratedTree(tree, times, [])
    order = {n.node_id: i for i, n in enumerate(solver.internal_nodes)}
    root_id = tree.node_id
    child_id = next(n.node_id for n in tree.walk() if n.parent and n.parent.node_id == root_id and not n.is_tip())

    baseline = np.array(
        [times.times[solver.internal_nodes[i].node_id] for i in range(len(order))]
    )
    solver._check_time_order(baseline)  # a sane vector must not raise

    negative = baseline.copy()
    negative[order[child_id]] = -1.0
    with pytest.raises(ValueError, match="older than the present") as exc:
        solver._check_time_order(negative)
    assert str(child_id) in str(exc.value)

    inverted = baseline.copy()
    inverted[order[child_id]] = baseline[order[root_id]] * 2.0
    with pytest.raises(ValueError, match="parent-child inversion") as exc:
        solver._check_time_order(inverted)
    assert f"node {child_id}" in str(exc.value)
    assert f"ancestor node {root_id}" in str(exc.value)


# ---------------------------------------------------------------------------
# G7: one-sided bound sets take their finite edge (unit-invariantly)
# ---------------------------------------------------------------------------

# A|B has relative time 0.23334 and the root 0.376 here, so a one-sided bound
# can sit either side of the uncalibrated scale and exercise both branches.
_TWO_SIDED = "((A:0.20,B:0.30):0.15,(C:0.10,D:0.25):0.20):0.0;"


def _two_sided_tree():
    return parse_newick(_TWO_SIDED)


@pytest.mark.gate("G7")
def test_one_sided_min_bound_below_the_uncalibrated_age_is_not_inert() -> None:
    """A minimum bound younger than the RRF solution pins the node to it.

    The factor used to be floored at ``1.0``, i.e. at the *uncalibrated* root
    age of this arbitrary relative scale, so this bound did nothing at all and
    the node stayed at 0.23334 with no warning.
    """
    tree = _two_sided_tree()
    times = rrf_times(tree)
    ab = tree.mrca(["A", "B"])
    cal = calibrate(times, [Calibration(node_id=ab.node_id, min_bound=0.1)], method="bounds")
    assert cal.times[ab.node_id] == pytest.approx(0.1, rel=1e-9)
    assert cal.warnings == []


@pytest.mark.gate("G7")
def test_one_sided_max_bound_above_the_uncalibrated_age_does_not_inflate() -> None:
    """A maximum bound older than the RRF solution pins downward, not upward."""
    tree = _two_sided_tree()
    times = rrf_times(tree)
    ab = tree.mrca(["A", "B"])
    uncalibrated = times.times[ab.node_id]
    cal = calibrate(times, [Calibration(node_id=ab.node_id, max_bound=uncalibrated / 2.0)],
                    method="bounds")
    assert cal.times[ab.node_id] == pytest.approx(uncalibrated / 2.0, rel=1e-9)
    # every age moved *with* the bound: the factor is a scale, not a floor
    assert cal.time_factor == pytest.approx(0.5, rel=1e-9)


@pytest.mark.gate("G7")
def test_one_sided_bounds_are_invariant_to_the_time_unit() -> None:
    """The same physical constraint in Ma and in Ga gives ages 1000x apart.

    Flooring the factor at the dimensionless 1.0 compared a time-unit-dependent
    quantity against a unit-free constant, so ``min 0.1`` (Ma) and ``min 100``
    (the same bound in Ga) took different code paths on the same tree.
    """
    ma = []
    for bound in (0.1, 100.0):
        tree = _two_sided_tree()
        cal = calibrate(
            rrf_times(tree),
            [Calibration(node_id=tree.mrca(["A", "B"]).node_id, min_bound=bound)],
            method="bounds",
        )
        ma.append(cal)
    ratio = [b / a for a, b in zip(ma[0].times.values(), ma[1].times.values()) if a > 0]
    assert all(r == pytest.approx(1000.0, rel=1e-9) for r in ratio), ratio
    assert ma[1].time_factor == pytest.approx(1000.0 * ma[0].time_factor, rel=1e-9)


# ---------------------------------------------------------------------------
# G7: bound-adjacent nodes must not be mislabelled, and off-scale input warned
# ---------------------------------------------------------------------------

@pytest.mark.gate("G7")
def test_node_on_its_bound_is_not_reported_as_unsatisfied() -> None:
    """A one-ULP distance to the bound is 'on the bound', not a failed solve.

    The truncation test and the pinned-node test next to it already compared
    with ``_PIN_TOL``; the 'constraint not satisfied' flag used an exact
    comparison, so a node the solver *had* placed on its maximum bound was
    written into the user-facing ``notes`` column as a constraint the solver
    could not reach.
    """
    tree = _simple_tree()
    times = rrf_times(tree, normalize=True)
    ab = tree.mrca(["A", "B"])
    cd = tree.mrca(["D", "E"])
    cal = calibrate(
        times,
        [
            Calibration(node_id=ab.node_id, min_bound=1.0, max_bound=2.0),
            Calibration(
                node_id=cd.node_id,
                min_bound=0.35,
                max_bound=0.55,
                density="exponential",
                density_params={"offset": 0.3, "mean": 0.5},
            ),
        ],
        method="bounds",
    )
    assert cal.times[cd.node_id] == pytest.approx(0.55)
    ci = confidence_interval(cal, seq_length=500)
    notes = ci.table.loc[ci.table["node_id"] == cd.node_id, "notes"].iloc[0]
    assert "constraint not satisfied" not in notes
    assert "pinned at its maximum hard bound" in notes
    # ...and the density quantiles the pinned interval came from contain it
    lo, hi = density_quantiles("exponential", {"offset": 0.3, "mean": 0.5})
    row = ci.table.loc[ci.table["node_id"] == cd.node_id].iloc[0]
    assert lo <= float(row["lower"]) <= float(row["upper"]) <= hi


@pytest.mark.gate("G7")
def test_genuinely_unsatisfied_constraint_is_still_flagged() -> None:
    """The tolerance must not silence the case the flag is written for."""
    tree = _simple_tree()
    times = rrf_times(tree, normalize=True)
    ab = tree.mrca(["A", "B"])
    cal = calibrate(
        times, [Calibration(node_id=ab.node_id, min_bound=1.0, max_bound=2.0)],
        method="bounds",
    )
    # the point estimate stays where the bounds it was solved with put it,
    # while the *reported* minimum is stretched far past it: the guard has to
    # keep saying so, because the printed interval is no eq-(7) statement
    unreachable = replace(
        cal,
        params={**cal.params, "method": "effective"},
        calibrations=[
            replace(
                cal.calibrations[0],
                min_bound=4.0,
                used_min_bound=1.0,
                used_max_bound=2.0,
            )
        ],
    )
    ci = confidence_interval(unreachable, seq_length=500)
    notes = ci.table.loc[ci.table["node_id"] == ab.node_id, "notes"].iloc[0]
    assert "constraint not satisfied" in notes


@pytest.mark.gate("G7")
def test_off_scale_branch_lengths_are_warned_once(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Tiny branch lengths widen the interval without the estimator changing.

    msz236 eq (14) is scale-free in the *ratio* but not in the resulting
    half-widths, so a tree not in substitutions per site silently gets a
    differently-sized interval.  The estimator is the published one and is left
    alone; the user is told instead.
    """
    tree = _simple_tree()
    for node in tree.walk():
        if node.blen is not None:
            node.blen *= 1e-4
    times = rrf_times(tree)
    ab = tree.mrca(["A", "B"])
    cal = calibrate(
        times, [Calibration(node_id=ab.node_id, min_bound=1.0, max_bound=2.0)],
        method="bounds",
    )
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="openreltime"):
        confidence_interval(cal, seq_length=1000)
    messages = [r.getMessage() for r in caplog.records if "substitutions" in r.getMessage()]
    assert len(messages) == 1, messages
    assert "longest branch length" in messages[0]

    # an in-scale tree says nothing
    ok = _simple_tree()
    ok_times = rrf_times(ok)
    ok_cal = calibrate(
        ok_times,
        [Calibration(node_id=ok.mrca(["A", "B"]).node_id, min_bound=1.0, max_bound=2.0)],
        method="bounds",
    )
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="openreltime"):
        confidence_interval(ok_cal, seq_length=1000)
    assert not [r for r in caplog.records if "substitutions" in r.getMessage()]
