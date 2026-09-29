"""End-to-end tests for CorrTest and ddBD."""

from __future__ import annotations

import math

import pytest

from e2e_tests.conftest import (  # noqa: E402
    GOLDEN_DIR,
    OUTGROUP_STR,
    TREE,
)


def test_corrtest_cli_matches_golden(run_cli, tmp_out):
    """CorrTest CLI reproduces the R3F golden score."""
    run_cli(
        "corrtest",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "-o", "ct",
    )
    txt = tmp_out / "ct_corrtest.txt"
    assert txt.exists()
    lines = txt.read_text().splitlines()
    score = float(next(ln for ln in lines if ln.startswith("score")).split("=")[1])
    p = next(ln for ln in lines if ln.startswith("P"))
    golden_score = float(
        (GOLDEN_DIR / "example_sr0_corrtest.txt").read_text()
        .splitlines()[0].split("=")[1]
    )
    assert score == pytest.approx(golden_score, abs=1e-4)
    assert p == "P-value < 0.001"


def test_corrtest_api_full_precision(api_tree):
    """API CorrScore equals the unrounded R3F double."""
    import openreltime as ort

    res = ort.corrtest(api_tree)
    exact = float(
        next(
            ln for ln in (GOLDEN_DIR / "example_corrtest_full_precision.txt").read_text().splitlines()
            if ln.startswith("sr0")
        ).split("\t")[1]
    )
    assert res.score == exact


def test_corrtest_resample(run_cli, tmp_out):
    """Sister resampling yields a reproducible score with a fixed seed."""
    run_cli(
        "corrtest",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "--sister-resample", "50",
        "--seed", "1",
        "-o", "ct_sr",
    )
    txt = tmp_out / "ct_sr_corrtest.txt"
    score = float(txt.read_text().splitlines()[0].split("=")[1])
    golden_sr50 = float(
        (GOLDEN_DIR / "example_sr50_corrtest.txt").read_text()
        .splitlines()[0].split("=")[1]
    )
    assert score == pytest.approx(golden_sr50, abs=1e-4)


def test_ddbd_cli_matches_golden(run_cli, tmp_out):
    """ddBD CLI reproduces the R3F no-anchor parameters."""
    run_cli(
        "ddbd",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "-o", "bd",
    )
    txt = tmp_out / "bd_ddbd.txt"
    assert txt.exists()
    lines = txt.read_text().splitlines()
    gb, gd, gs = (float(v) for v in lines[1].split("\t"))
    golden = (
        GOLDEN_DIR / "example_noanchor_ddbd.txt"
    ).read_text().splitlines()[1].split("\t")
    assert gb == pytest.approx(float(golden[0]), rel=5e-4)
    assert gd == pytest.approx(float(golden[1]), rel=5e-4)
    assert gs == pytest.approx(float(golden[2]), rel=5e-3)


def test_ddbd_anchor_matches_golden(run_cli, tmp_out):
    """Anchored ddBD reproduces the R3F anchor run."""
    run_cli(
        "ddbd",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "--anchor-node", "272",
        "--anchor-time", "1.85",
        "-o", "bd_anchor",
    )
    txt = tmp_out / "bd_anchor_ddbd.txt"
    lines = txt.read_text().splitlines()
    gb, gd, gs = (float(v) for v in lines[1].split("\t"))
    golden = (
        GOLDEN_DIR / "example_anchor_ddbd.txt"
    ).read_text().splitlines()[1].split("\t")
    assert gb == pytest.approx(float(golden[0]), rel=5e-4)
    assert gd == pytest.approx(float(golden[1]), rel=5e-4)
    assert gs == pytest.approx(float(golden[2]), rel=5e-3)


def test_ddbd_plot(run_cli, tmp_out):
    """ddBD --plot writes a non-empty PNG."""
    plot = tmp_out / "bd_curve.png"
    run_cli(
        "ddbd",
        "-i", str(TREE),
        "--outgroup", OUTGROUP_STR,
        "--plot", str(plot),
        "-o", "bd_plot",
    )
    assert plot.exists() and plot.stat().st_size > 0


def test_ddbd_api_matches_cli(api_tree, tmp_out):
    """API ddBD returns the same birth/death/sampling as a CLI run."""
    import openreltime as ort

    api_res = ort.ddbd(api_tree)
    assert math.isfinite(api_res.birth_rate)
    assert math.isfinite(api_res.death_rate)
    assert 0 < api_res.sampling_frac <= 1
