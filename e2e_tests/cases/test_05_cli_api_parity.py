"""End-to-end tests for CLI/API parity and report.json merging."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from e2e_tests.conftest import (  # noqa: E402
    OUTGROUP_STR,
    TREE,
)

HARD_BOUNDS_CALIBRATIONS = (
    "node_id\ttaxon_set\tmin_bound\tmax_bound\tdensity\tdensity_params\n"
    "-\tHomo_sapiens|Pan_troglodytes\t6.0\t8.5\t.\t.\n"
    "-\tDasypus_novemcinctus|Choloepus_didactylus\t70.0\t105.0\t.\t.\n"
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


def _api_write(sub: str, prefix: Path, calibrations: Path) -> None:
    """Run *sub* through the Python API and write its artefacts at *prefix*.

    The flags mirror the CLI subcommand's own ``write()`` call exactly, so any
    difference in the emitted bytes is a real divergence between the two entry
    points rather than a difference in what they were asked to write.
    """
    import openreltime as ort

    tree = ort.read_tree(str(TREE), outgroup=OUTGROUP_STR.split(","))
    if sub == "rates":
        ort.rrf_rates(tree).write(prefix)
    elif sub == "times":
        ort.rrf_times(tree).write(prefix, with_rate=False)
    elif sub == "rates-times":
        ort.rrf_rates_times(tree).write(prefix, with_rate=True, nexus=True)
    elif sub == "calibrate":
        cals = ort.parse_calibrations(str(calibrations))
        ort.calibrate(ort.rrf_times(tree), cals).write(prefix)
    else:  # pragma: no cover - the parameter list is closed
        raise AssertionError(sub)


@pytest.mark.parametrize(
    "sub", ["rates", "times", "rates-times", "calibrate"]
)
def test_cli_and_api_artefacts_are_byte_identical_on_one_machine(
    run_cli, tmp_out, sub
):
    """Gate G3, byte level: the API and the CLI emit the same bytes.

    The numeric comparison above (`np.allclose` at 1e-9) only checks two
    columns of one CSV.  Both entry points call the same ``report.write()``, so
    every artefact a subcommand writes -- the node tables, the annotated
    timetrees, the calibrated and effective-bound CSVs -- must agree byte for
    byte, and a divergence in column order, float representation or Newick
    precision must fail here rather than surface as two irreconcilable results
    files.

    Scope is deliberately narrow: same machine, same interpreter, one process
    for the API and one for the CLI, deterministic options (no ``--resolve
    random``, so no seed is needed).  It is *not* a cross-platform claim -- CSV
    line endings, locale and float reprs may differ between operating systems.

    The JSON report is compared as a document rather than as bytes: the CLI
    records tree-loading provenance (``resolve_polytomy``, ``seed``, the
    resolved input paths) that the API has no reason to invent, so the
    assertion is that the two blocks agree wherever they both state a value and
    that the CLI block adds only those extra keys.
    """
    calibrations = tmp_out / "hard_bounds.tsv"
    calibrations.write_text(HARD_BOUNDS_CALIBRATIONS)
    api_prefix = tmp_out / "api"
    cli_prefix = tmp_out / "cli"

    _api_write(sub, api_prefix, calibrations)
    cli_args = [sub, "-i", str(TREE), "--outgroup", OUTGROUP_STR, "-o", str(cli_prefix)]
    if sub == "calibrate":
        cli_args[3:3] = ["-c", str(calibrations)]
    run_cli(*cli_args)

    api_files = {p.name.removeprefix("api_"): p for p in tmp_out.glob("api_*")}
    cli_files = {p.name.removeprefix("cli_"): p for p in tmp_out.glob("cli_*")}
    assert set(api_files) == set(cli_files), (
        f"{sub}: the two entry points wrote different artefact sets"
    )
    for name in sorted(api_files):
        if name == "report.json":
            continue
        api_bytes = api_files[name].read_bytes()
        cli_bytes = cli_files[name].read_bytes()
        assert api_bytes == cli_bytes, f"{sub}: {name} differs byte-for-byte"

    api_report = json.loads(api_files["report.json"].read_text())
    cli_report = json.loads(cli_files["report.json"].read_text())
    assert set(api_report) == set(cli_report)
    for block in api_report:
        api_params = api_report[block]["parameters"]
        cli_params = cli_report[block]["parameters"]
        shared = set(api_params) & set(cli_params)
        assert shared, f"{block}: the two reports share no parameter"
        for key in shared:
            assert api_params[key] == cli_params[key], f"{block}.{key} differs"
        assert api_report[block].get("warnings") == cli_report[block].get("warnings")


def test_report_json_merges_blocks(run_cli, tmp_out, effective_calibrations):
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
        "-c", str(effective_calibrations),
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
