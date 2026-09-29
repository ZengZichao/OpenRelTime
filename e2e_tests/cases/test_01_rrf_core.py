"""End-to-end tests for the RRF core commands (rates, times, rates-times)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))
from compare_golden import compare  # noqa: E402

from e2e_tests.conftest import (  # noqa: E402
    GOLDEN_DIR,
    OUTGROUP,
    OUTGROUP_STR,
    TREE,
)


def test_rates_cli_regression(run_cli, tmp_out, read_report):
    """`rates` CLI reproduces R3F golden rates."""
    run_cli(
        "rates",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "-o", "rates",
    )
    csv = tmp_out / "rates_rates.csv"
    nwk = tmp_out / "rates_rates.nwk"
    assert csv.exists() and csv.stat().st_size > 0
    assert nwk.exists() and nwk.stat().st_size > 0

    result = compare(TREE, GOLDEN_DIR, OUTGROUP, with_rate=True)
    assert result["rate_slope"] >= 0.999
    assert result["rate_median_rel"] <= 1e-6

    report = read_report("rates", analysis="rrf_rates")
    assert report["analysis"] == "rrf_rates"
    assert report["parameters"]["mean"] == "geometric"


def test_times_cli_r3f_compat(run_cli, tmp_out, read_report):
    """``times`` CLI writes R3F-compatible files and reproduces golden times."""
    run_cli(
        "times",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "--r3f-compat",
        "-o", "times",
    )
    expected = {
        tmp_out / "times_RRF_times.csv",
        tmp_out / "times_RRF_timetree.nwk",
        tmp_out / "times_report.json",
    }
    for path in expected:
        assert path.exists() and path.stat().st_size > 0, f"missing {path.name}"

    result = compare(TREE, GOLDEN_DIR, OUTGROUP, with_rate=False)
    assert result["time_slope"] >= 0.999
    assert result["time_max_rel"] <= 1e-6
    assert result["n_common"] == 541

    report = read_report("times", analysis="rrf_times")
    assert report["analysis"] == "rrf_times"
    # r3f_compat affects file naming but is not recorded in the report
    assert (tmp_out / "times_RRF_times.csv").exists()


def test_rates_times_cli_outputs(run_cli, tmp_out, read_report):
    """``rates-times`` produces the full file set and a valid NEXUS."""
    plot_path = tmp_out / "rt_timetree.png"
    run_cli(
        "rates-times",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "--plot", str(plot_path),
        "-o", "rt",
    )
    expected = {
        tmp_out / "rt_rates_times.csv",
        tmp_out / "rt_timetree.nexus",
        tmp_out / "rt_report.json",
        plot_path,
    }
    for path in expected:
        assert path.exists() and path.stat().st_size > 0, f"missing {path.name}"

    nexus = (tmp_out / "rt_timetree.nexus").read_text()
    assert "[&rate=" in nexus

    report = read_report("rt", analysis="rrf_times")
    assert report["analysis"] == "rrf_times"


def test_times_normalize(run_cli, tmp_out):
    """``--normalize`` scales the root age to 1.0."""
    run_cli(
        "times",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "--normalize",
        "-o", "norm",
    )
    import pandas as pd
    df = pd.read_csv(tmp_out / "norm_times.csv")
    assert df["Time"].max() == pytest.approx(1.0, abs=1e-9)


def test_arithmetic_mean_cli(run_cli, tmp_out, read_report):
    """``--mean arithmetic`` is accepted and recorded."""
    run_cli(
        "times",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "--mean", "arithmetic",
        "-o", "arith",
    )
    report = read_report("arith", analysis="rrf_times")
    assert report["parameters"]["mean"] == "arithmetic"


def test_guard_options(run_cli, tmp_out, read_report, api_tree):
    """Guard threshold/no-guard options reach the report."""
    run_cli(
        "times",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "--rate-ratio-threshold", "10",
        "-o", "guard10",
    )
    report = read_report("guard10", analysis="rrf_times")
    assert report["parameters"]["rate_ratio_threshold"] == 10.0

    run_cli(
        "times",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "--no-guard",
        "-o", "noguard",
    )
    report = read_report("noguard", analysis="rrf_times")
    assert report["parameters"]["rate_ratio_threshold"] is None

    # API parity for no-guard
    import openreltime as ort
    res = ort.rrf_times(api_tree, rate_ratio_threshold=None)
    assert all(t >= 0 and t == t for t in res.times.values() if t is not None)


def test_nexus_roundtrip(run_cli, tmp_out):
    """Rates annotated in NEXUS survive a read/write roundtrip."""
    run_cli(
        "rates-times",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "-o", "rt_for_nexus",
    )
    import openreltime as ort

    nexus_in = tmp_out / "rt_for_nexus_timetree.nexus"
    t = ort.read_tree(str(nexus_in), fmt="nexus")
    rate_map = {n.node_id: n.rate for n in t.walk() if hasattr(n, "rate")}
    out_nexus = tmp_out / "roundtrip.nexus"
    nexus_text = ort.write_nexus(t, rate_map=rate_map)
    out_nexus.write_text(nexus_text)
    assert "[&rate=" in nexus_text
