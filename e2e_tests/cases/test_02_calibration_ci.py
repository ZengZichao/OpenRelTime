"""End-to-end tests for calibration and analytical confidence intervals."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from e2e_tests.conftest import (  # noqa: E402
    CALIBRATIONS,
    OUTGROUP_STR,
    TREE,
)


@pytest.fixture
def hard_bounds_cals(tmp_out: Path) -> Path:
    """A calibration file with only hard bounds (no densities)."""
    path = tmp_out / "hard_bounds.tsv"
    path.write_text(
        "node_id\ttaxon_set\tmin_bound\tmax_bound\tdensity\tdensity_params\n"
        "-\tHomo_sapiens|Pan_troglodytes\t6.0\t8.5\t.\t.\n"
        "-\tDasypus_novemcinctus|Choloepus_didactylus\t70.0\t105.0\t.\t.\n"
    )
    return path


def test_calibrate_bounds_cli(run_cli, tmp_out, read_report, hard_bounds_cals):
    """Hard-bound calibration produces absolute ages and a report."""
    run_cli(
        "calibrate",
        "-i", str(TREE),
        "-c", str(hard_bounds_cals),
        "--outgroup", OUTGROUP_STR,
        "--method", "bounds",
        "-o", "cal_bounds",
    )
    expected = {
        tmp_out / "cal_bounds_calibrated.csv",
        tmp_out / "cal_bounds_calibrated.nwk",
        tmp_out / "cal_bounds_calibrated.nexus",
        tmp_out / "cal_bounds_report.json",
    }
    for path in expected:
        assert path.exists() and path.stat().st_size > 0

    df = pd.read_csv(tmp_out / "cal_bounds_calibrated.csv")
    assert (df["Time"] >= 0).all()
    assert df["Time"].max() > 0

    report = read_report("cal_bounds", analysis="calibrate")
    assert report["parameters"]["method"] == "bounds"
    assert report["parameters"]["mean"] == "geometric"


def test_calibrate_effective_cli(run_cli, tmp_out, read_report):
    """Effective-bounds calibration samples densities and writes extra output."""
    run_cli(
        "calibrate",
        "-i", str(TREE),
        "-c", str(CALIBRATIONS),
        "--outgroup", OUTGROUP_STR,
        "--method", "effective",
        "--n-effective", "5000",
        "--seed", "42",
        "-o", "cal_eff",
    )
    expected = {
        tmp_out / "cal_eff_calibrated.csv",
        tmp_out / "cal_eff_effective_bounds.csv",
        tmp_out / "cal_eff_report.json",
    }
    for path in expected:
        assert path.exists() and path.stat().st_size > 0

    report = read_report("cal_eff", analysis="calibrate")
    assert report["parameters"]["method"] == "effective"
    assert report["parameters"]["n_effective"] == 5000
    assert report["parameters"]["seed"] == 42


def test_ci_default(run_cli, tmp_out, read_report, hard_bounds_cals):
    """CI from a previous calibration with rate-heterogeneity variance only."""
    run_cli(
        "calibrate",
        "-i", str(TREE),
        "-c", str(hard_bounds_cals),
        "--outgroup", OUTGROUP_STR,
        "-o", "cal_ci",
    )
    run_cli(
        "ci",
        "-c", "cal_ci",
        "-o", "ci_default",
    )
    ci_csv = tmp_out / "ci_default_ci.csv"
    assert ci_csv.exists()
    df = pd.read_csv(ci_csv)
    assert (df["time"] >= 0).all()
    assert (df["lower"] <= df["time"]).all()
    assert (df["time"] <= df["upper"]).all()
    assert (df["upper"] - df["lower"] >= 0).all()

    report = read_report("ci_default", analysis="confidence_interval")
    assert report["level"] == 0.95
    assert "method=bounds" in report["method"]


def test_ci_with_n_sites(run_cli, tmp_out, read_report, hard_bounds_cals):
    """CI widens when a small sequence length is supplied."""
    run_cli(
        "calibrate",
        "-i", str(TREE),
        "-c", str(hard_bounds_cals),
        "--outgroup", OUTGROUP_STR,
        "-o", "cal_ns",
    )
    run_cli("ci", "-c", "cal_ns", "-o", "ci_vh")
    run_cli("ci", "-c", "cal_ns", "--n-sites", "10", "-o", "ci_ns")

    df_vh = pd.read_csv(tmp_out / "ci_vh_ci.csv")
    df_ns = pd.read_csv(tmp_out / "ci_ns_ci.csv")
    # align by node_id
    merged = pd.merge(df_vh, df_ns, on="node_id", suffixes=("_vh", "_ns"))
    wider = (merged["width_ns"] >= merged["width_vh"] - 1e-12).mean()
    assert wider >= 0.95  # most nodes should be wider with tiny L

    report = read_report("ci_ns", analysis="confidence_interval")
    assert "L=10" in report["method"]


def test_ci_with_branch_var(run_cli, tmp_out, hard_bounds_cals):
    """CI responds to a user-supplied per-edge sampling variance."""
    run_cli(
        "calibrate",
        "-i", str(TREE),
        "-c", str(hard_bounds_cals),
        "--outgroup", OUTGROUP_STR,
        "-o", "cal_bv",
    )
    # build branch-var with tiny variance for all internal edges
    cal = pd.read_csv(tmp_out / "cal_bv_calibrated.csv")
    internal = cal[cal["NodeId"] > len(cal) // 2].copy()
    internal = internal[["NodeId"]].rename(columns={"NodeId": "node_id"})
    internal["var"] = 1e-6
    bv_path = tmp_out / "branch_var.tsv"
    internal.to_csv(bv_path, sep="\t", index=False)

    run_cli("ci", "-c", "cal_bv", "--branch-var", str(bv_path), "-o", "ci_bv")
    df = pd.read_csv(tmp_out / "ci_bv_ci.csv")
    assert (df["lower"] <= df["time"]).all()
    assert (df["time"] <= df["upper"]).all()


def test_calibration_api_matches_cli(run_cli, tmp_out, api_tree, hard_bounds_cals):
    """Python API calibration matches the CLI output for the same settings."""
    import openreltime as ort

    cals = ort.parse_calibrations(str(hard_bounds_cals))
    times = ort.rrf_times(api_tree)
    api_result = ort.calibrate(times, cals, method="bounds")
    api_result.write(str(tmp_out / "api_cal"))

    run_cli(
        "calibrate",
        "-i", str(TREE),
        "-c", str(hard_bounds_cals),
        "--outgroup", OUTGROUP_STR,
        "--method", "bounds",
        "-o", "cli_cal",
    )
    api_df = pd.read_csv(tmp_out / "api_cal_calibrated.csv")
    cli_df = pd.read_csv(tmp_out / "cli_cal_calibrated.csv")
    api_df = api_df.sort_values("NodeId").reset_index(drop=True)
    cli_df = cli_df.sort_values("NodeId").reset_index(drop=True)
    assert np.allclose(api_df["Time"], cli_df["Time"], rtol=1e-9, atol=1e-9)


def test_conflicting_calibration_rejected(run_cli_raw, tmp_out):
    """A calibration row with min > max raises a clean error."""
    bad = tmp_out / "bad_cals.tsv"
    bad.write_text(
        "node_id\ttaxon_set\tmin_bound\tmax_bound\tdensity\tdensity_params\n"
        "-\tHomo_sapiens|Pan_troglodytes\t8.5\t6.0\t.\t.\n"
    )
    proc = run_cli_raw(
        "calibrate",
        "-i", str(TREE),
        "-c", str(bad),
        "--outgroup", OUTGROUP_STR,
        "-o", "bad_cal",
    )
    assert proc.returncode != 0
    assert "min" in proc.stderr.lower() and "max" in proc.stderr.lower()


def test_density_only_requires_effective(run_cli_raw, tmp_out):
    """Density-only row with --method bounds is rejected."""
    bad = tmp_out / "density_only.tsv"
    bad.write_text(
        "node_id\ttaxon_set\tmin_bound\tmax_bound\tdensity\tdensity_params\n"
        "-\tHomo_sapiens|Pan_troglodytes\t.\t.\texponential\toffset=60;mean=20\n"
    )
    proc = run_cli_raw(
        "calibrate",
        "-i", str(TREE),
        "-c", str(bad),
        "--outgroup", OUTGROUP_STR,
        "--method", "bounds",
        "-o", "density_bounds",
    )
    assert proc.returncode != 0
    assert "effective" in proc.stderr.lower()


def test_pinned_node_ci_flagged(run_cli, tmp_out):
    """A node pinned at its own hard bound is flagged se_reliable=False."""
    # The MRCA of Homo and Pan calibrates near 7.25 under the unbounded RRF times;
    # pin it with a min/max bracket that is much tighter than the delta-method SE.
    tight = tmp_out / "tight.tsv"
    tight.write_text(
        "node_id\ttaxon_set\tmin_bound\tmax_bound\tdensity\tdensity_params\n"
        "-\tHomo_sapiens|Pan_troglodytes\t7.25\t7.25001\t.\t.\n"
    )
    run_cli(
        "calibrate",
        "-i", str(TREE),
        "-c", str(tight),
        "--outgroup", OUTGROUP_STR,
        "-o", "pin_cal",
    )
    run_cli("ci", "-c", "pin_cal", "-o", "pin_ci")
    df = pd.read_csv(tmp_out / "pin_ci_ci.csv")
    row = df[df["node_id"] == 371]
    assert not row.empty
    assert row.iloc[0]["se_reliable"] == False
    assert "calibration" in row.iloc[0]["notes"].lower()
