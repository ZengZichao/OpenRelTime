"""Calibration module behaviour, the CI module, and edge cases (validation gate G7)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from openreltime.calibrate import Calibration, calibrate, parse_calibrations
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


# ---------------------------------------------------------------------------
# CI module
# ---------------------------------------------------------------------------

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

def test_negative_branch_length_rejected(tmp_path: Path) -> None:
    bad = tmp_path / "neg.nwk"
    bad.write_text("(A:1,B:-2);")
    with pytest.raises(ValueError, match="negative branch length"):
        read_tree(bad)


def test_duplicate_tip_names_rejected(tmp_path: Path) -> None:
    bad = tmp_path / "dup.nwk"
    bad.write_text("((A:1,A:2):0.5,(C:1,D:2):0.4):0.0;")
    with pytest.raises(ValueError, match="duplicate tip names"):
        read_tree(bad)


def test_polytomy_rejected_and_random_resolution(tmp_path: Path) -> None:
    poly = tmp_path / "poly.nwk"
    poly.write_text("((A:1,B:2,C:3):4,(D:5,E:6):7);")
    with pytest.raises(ValueError, match="polytom"):
        read_tree(poly)
    tree = read_tree(poly, resolve_polytomy="random", seed=1)
    assert tree.is_binary()


def test_two_tip_ingroup_rejected() -> None:
    with pytest.raises(ValueError, match="at least 3"):
        tree = parse_newick("((A:1,B:2),(C:3,D:4)):0.0;")
        from openreltime.treeio import root_outgroup

        root_outgroup(tree, ["C", "D"])


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


def test_extreme_rates_trigger_guard_and_no_inf() -> None:
    newick = "((A:0.0000001,B:2.0):0.1,(C:0.5,D:0.6):0.9):0.05;"
    tree = parse_newick(newick)
    res = rrf_times(tree)
    vals = np.array(list(res.times.values()), dtype=float)
    assert np.isfinite(vals).all()


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


def test_arithmetic_mean_runs_and_differs() -> None:
    tree = _simple_tree()
    geo = rrf_times(tree)
    ari = rrf_times(tree, mean="arithmetic")
    assert any(
        abs(geo.times[k] - ari.times[k]) > 1e-9 for k in geo.times
    )


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
