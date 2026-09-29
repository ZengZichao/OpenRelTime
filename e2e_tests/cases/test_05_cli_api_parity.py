"""End-to-end tests for CLI/API parity and report.json merging."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from e2e_tests.conftest import (  # noqa: E402
    CALIBRATIONS,
    OUTGROUP_STR,
    TREE,
)


def test_api_write_matches_cli_rates_times(run_cli, tmp_out):
    """API ``write()`` produces the same files and numbers as the CLI."""
    import openreltime as ort

    tree = ort.read_tree(str(TREE), outgroup=OUTGROUP_STR.split(","))
    api_res = ort.rrf_rates_times(tree)
    api_written = api_res.write(str(tmp_out / "api_rt"), with_rate=True, nexus=True)

    run_cli(
        "rates-times",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "-o", "cli_rt",
    )
    api_csv = next(p for p in api_written if p.suffix == ".csv")
    cli_csv = tmp_out / "cli_rt_rates_times.csv"

    api_df = pd.read_csv(api_csv).sort_values("NodeId").reset_index(drop=True)
    cli_df = pd.read_csv(cli_csv).sort_values("NodeId").reset_index(drop=True)
    assert np.allclose(api_df["Time"], cli_df["Time"], rtol=1e-9, atol=1e-9)
    assert np.allclose(api_df["Rate"], cli_df["Rate"], rtol=1e-9, atol=1e-9)


def test_report_json_merges_blocks(run_cli, tmp_out):
    """Multiple analyses with the same prefix merge into one report."""
    run_cli(
        "rates",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "-o", "merged",
    )
    run_cli(
        "times",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "-o", "merged",
    )
    run_cli(
        "calibrate",
        "-i", str(TREE),
        "-c", str(CALIBRATIONS),
        "--outgroup", OUTGROUP_STR,
        "--method", "effective",
        "--n-effective", "1000",
        "--seed", "42",
        "-o", "merged",
    )
    run_cli(
        "ci",
        "-c", "merged",
        "-o", "merged",
    )
    report_path = tmp_out / "merged_report.json"
    assert report_path.exists()
    data = json.loads(report_path.read_text())
    assert set(data.keys()) == {"rrf_rates", "rrf_times", "calibrate", "confidence_interval"}
    assert data["rrf_rates"]["parameters"]["mean"] == "geometric"
    assert data["rrf_times"]["parameters"]["mean"] == "geometric"
    assert data["calibrate"]["parameters"]["method"] == "effective"
    assert data["confidence_interval"]["level"] == 0.95


def test_cli_parameters_recorded(run_cli, tmp_out, read_report):
    """Non-default CLI options are recorded in the report."""
    run_cli(
        "times",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "--mean", "arithmetic",
        "--no-guard",
        "--normalize",
        "-o", "params",
    )
    report = read_report("params", analysis="rrf_times")
    params = report["parameters"]
    assert params["mean"] == "arithmetic"
    assert params["rate_ratio_threshold"] is None
    assert params["normalize"] is True
