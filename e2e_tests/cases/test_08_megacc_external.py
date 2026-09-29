"""MEGA-CC bridge end-to-end tests (skipped unless megacc is available)."""

from __future__ import annotations

import os
import shutil

import pytest

from e2e_tests.conftest import (  # noqa: E402
    CALIBRATIONS,
    OUTGROUP_STR,
    TREE,
)


def _megacc_path() -> str | None:
    """Return the megacc executable, or ``None`` when it is not installed.

    ``None`` rather than ``pytest.skip``: see ``test_07_blb_external`` — a skip
    raised during import aborts collection for the whole suite.
    """
    env = os.environ.get("MEGACC")
    if env:
        return env
    return shutil.which("megacc")


pytestmark = pytest.mark.skipif(
    _megacc_path() is None,
    reason="MEGA-CC megacc not available; set MEGACC env var to enable",
)


def test_megacc_with_calibrations(run_cli, tmp_out):
    """megacc bridge runs and parses MEGA-CC output."""
    run_cli(
        "megacc",
        "-i", str(TREE),
        "-c", str(CALIBRATIONS),
        "--megacc", _megacc_path(),
        "--outgroup", OUTGROUP_STR,
        "-o", "mega",
    )
    csv = tmp_out / "mega_megacc_times.csv"
    assert csv.exists()

    import pandas as pd
    df = pd.read_csv(csv)
    assert "source" in df.columns
    assert not df.empty


def test_megacc_without_calibrations(run_cli, tmp_out):
    """megacc can run branch-length-only RelTime."""
    run_cli(
        "megacc",
        "-i", str(TREE),
        "--megacc", _megacc_path(),
        "--outgroup", OUTGROUP_STR,
        "-o", "mega_nc",
    )
    csv = tmp_out / "mega_nc_megacc_times.csv"
    assert csv.exists()
