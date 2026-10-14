"""Shared fixtures for the OpenRelTime end-to-end test suite."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

# Project paths
E2E_ROOT = Path(__file__).resolve().parent
DATA_DIR = E2E_ROOT / "data"
GOLDEN_DIR = DATA_DIR / "golden"
OUTPUT_DIR = E2E_ROOT / "outputs"

TREE = DATA_DIR / "example.nwk"
CALIBRATIONS = DATA_DIR / "example_calibrations.tsv"
TAXA = DATA_DIR / "example_taxa.tsv"
# Small synthetic alignment used only for the BLB/IQ-TREE smoke test.
# It is intentionally tiny so the test finishes in seconds while still exercising
# the full IQ-TREE subprocess pipeline.
ALIGNMENT = DATA_DIR / "small.fasta"

OUTGROUP = [
    "Ornithorhynchus_anatinus",
    "Zaglossus_bruijni",
    "Tachyglossus_aculeatus",
]
OUTGROUP_STR = ",".join(OUTGROUP)


def _child_env() -> dict[str, str]:
    """Environment for CLI subprocesses, with the repo root importable.

    Both CLI fixtures launch with ``cwd=`` a scratch directory, so a relative
    ``PYTHONPATH=.`` from the invoking shell never reaches the child and every
    test dies with ``ModuleNotFoundError: openreltime``.  CI is unaffected
    because it installs the distribution first; a bare source checkout is not.
    Injecting the absolute root makes the suite runnable either way.
    """
    env = dict(os.environ)
    root = str(E2E_ROOT.parent)
    existing = env.get("PYTHONPATH") or ""
    parts = [p for p in existing.split(os.pathsep) if p]
    if root not in parts:
        parts.insert(0, root)
    env["PYTHONPATH"] = os.pathsep.join(parts)
    return env


def pytest_configure(config):
    """Ensure output directory exists once per session."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


@pytest.fixture(scope="session")
def project_root() -> Path:
    return E2E_ROOT.parent


@pytest.fixture(scope="session")
def data_dir() -> Path:
    return DATA_DIR


@pytest.fixture(scope="session")
def golden_dir() -> Path:
    return GOLDEN_DIR


@pytest.fixture
def tmp_out(tmp_path: Path) -> Path:
    """Per-test output directory under e2e_tests/outputs/."""
    out = OUTPUT_DIR / tmp_path.name
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    return out


@pytest.fixture
def run_cli(tmp_out: Path):
    """Run an ``openreltime`` CLI command and return CompletedProcess."""

    def _run(*args: str, check: bool = True, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
        cmd = [sys.executable, "-m", "openreltime.cli", *args]
        proc = subprocess.run(
            cmd,
            cwd=str(cwd or tmp_out),
            capture_output=True,
            text=True,
            env=_child_env(),
        )
        if check and proc.returncode != 0:
            raise AssertionError(
                f"CLI failed (exit {proc.returncode}): {' '.join(cmd)}\n"
                f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
            )
        return proc

    return _run


@pytest.fixture
def run_cli_raw(tmp_out: Path):
    """Run CLI without raising on non-zero exit; useful for error tests."""

    def _run(*args: str, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
        cmd = [sys.executable, "-m", "openreltime.cli", *args]
        return subprocess.run(
            cmd,
            cwd=str(cwd or tmp_out),
            capture_output=True,
            text=True,
            env=_child_env(),
        )

    return _run


@pytest.fixture
def read_report(tmp_out: Path):
    """Read the report JSON produced by a command prefix.

    Returns the analysis block identified by *analysis* when given;
    otherwise returns the whole report dict.
    """

    def _read(prefix: str, analysis: str | None = None) -> dict[str, Any]:
        path = tmp_out / f"{prefix}_report.json"
        if not path.exists():
            raise FileNotFoundError(f"report not found: {path}")
        data = json.loads(path.read_text())
        if analysis is None:
            return data
        if analysis not in data:
            raise KeyError(f"analysis {analysis!r} not in report; keys: {list(data)}")
        return data[analysis]

    return _read


#: The rows of :data:`EFFECTIVE_CALIBRATION_TEXT` -- a density calibration the
#: ``effective`` method can actually solve.
EFFECTIVE_CALIBRATION_TEXT = (
    "node_id\ttaxon_set\tmin_bound\tmax_bound\tdensity\tdensity_params\n"
    "-\tHomo_sapiens|Pan_troglodytes\t5.0\t8.0\t.\t.\n"
    "-\tLoxodonta_africana|Elephas_maximus\t.\t.\tnormal\tmean=2.7;sd=0.3\n"
    "-\tMus_musculus|Rattus_norvegicus\t12.0\t18.0\t.\t.\n"
    "-\tDasypus_novemcinctus|Choloepus_didactylus\t43.0\t52.0\t.\t.\n"
)


@pytest.fixture
def effective_calibrations(tmp_out: Path) -> Path:
    """A density calibration table ``method='effective'`` can solve.

    ``CALIBRATIONS`` (``data/example_calibrations.tsv``) exists to illustrate
    the file format, and its placeholder bounds are mutually inconsistent: the
    ``Elephas``/``Loxodonta`` density implies a much older
    ``Dasypus``/``Choloepus`` node than that row's own 70-105 bound allows, so
    under ``method='effective'`` the solver converges on ages that are not
    nested in the tree -- a node at -4,000 Ma and four parent-child inversions
    -- which :func:`openreltime.calibrate.calibrate` now rejects.  These four
    rows are the feasible counterpart: three hard bounds plus one narrow
    ``normal`` density centred on the age the other three rows already imply,
    so a sampled pair never empties the feasible intersection.
    """
    path = tmp_out / "effective_calibrations.tsv"
    path.write_text(EFFECTIVE_CALIBRATION_TEXT)
    return path


@pytest.fixture
def iqtree_path() -> str | None:
    """Path to IQ-TREE executable, if available."""
    env = os.environ.get("IQTREE")
    if env:
        return env
    for name in ("iqtree2", "iqtree"):
        p = shutil.which(name)
        if p:
            return p
    return None


@pytest.fixture
def megacc_path() -> str | None:
    """Path to MEGA-CC megacc executable, if available."""
    env = os.environ.get("MEGACC")
    if env:
        return env
    p = shutil.which("megacc")
    return p


@pytest.fixture(scope="session")
def api_tree():
    """Load the example tree once per session via the public API."""
    import openreltime as ort

    return ort.read_tree(str(TREE), outgroup=OUTGROUP)
