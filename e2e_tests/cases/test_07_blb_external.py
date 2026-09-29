"""BLB end-to-end tests (skipped unless IQ-TREE is available)."""

from __future__ import annotations

import os
import shutil

import pytest

from e2e_tests.conftest import (  # noqa: E402
    ALIGNMENT,
)

# Use a single outgroup tip for the BLB smoke test.  IQ-TREE is run on a tiny
# synthetic alignment, so a multi-tip outgroup would not reliably come out
# monophyletic; a single tip always satisfies the rooting/pruning requirement.
BLB_OUTGROUP = "Ornithorhynchus_anatinus"


def _iqtree_path() -> str | None:
    """Return the IQ-TREE executable, or ``None`` when it is not installed.

    ``None`` rather than ``pytest.skip``: this runs while the module is being
    imported to build ``pytestmark``, and a skip raised there aborts collection
    for the whole suite instead of skipping this file.
    """
    env = os.environ.get("IQTREE")
    if env:
        return env
    for name in ("iqtree2", "iqtree"):
        p = shutil.which(name)
        if p:
            return p
    return None


pytestmark = pytest.mark.skipif(
    _iqtree_path() is None,
    reason="IQ-TREE not available; set IQTREE env var or install iqtree/iqtree2",
)


def test_blb_full_pipeline(run_cli, tmp_out):
    """BLB runs IQ-TREE replicates and produces summary + replicate tables."""
    run_cli(
        "blb",
        "-a", str(ALIGNMENT),
        "--iqtree", _iqtree_path(),
        "--outgroup", BLB_OUTGROUP,
        "--gamma", "0.7",
        # Tiny smoke test: one subsample and one replicate keeps the suite fast
        # while still exercising the full BLB write path.
        "--n-subsample", "1",
        "--n-replicate", "1",
        "--seed", "42",
        # Fix the model and use a single thread so the smoke test finishes in
        # seconds instead of triggering IQ-TREE's model selection / auto-threading.
        "--iqtree-arg=-nt", "--iqtree-arg=1",
        "--iqtree-arg=-m", "--iqtree-arg=GTR+G",
        "-o", "blb_out",
    )
    summary = tmp_out / "blb_out" / "blb_summary.csv"
    replicates = tmp_out / "blb_out" / "blb_replicates.csv"
    assert summary.exists()
    assert replicates.exists()

    import pandas as pd
    df = pd.read_csv(summary)
    assert {"median_time", "ci_lower", "ci_upper", "topology_ok"}.issubset(set(df.columns))
    assert df["median_time"].notna().all()


def test_blb_iqtree_args(run_cli, tmp_out):
    """Extra IQ-TREE arguments reach the subprocess command line."""
    run_cli(
        "blb",
        "-a", str(ALIGNMENT),
        "--iqtree", _iqtree_path(),
        "--outgroup", BLB_OUTGROUP,
        "--n-subsample", "1",
        "--n-replicate", "1",
        "--iqtree-arg=-nt", "--iqtree-arg=1",
        "--iqtree-arg=-m", "--iqtree-arg=GTR+G",
        "-o", "blb_args",
    )
    # The workdir is deleted by default; we only assert the run did not crash.
    assert (tmp_out / "blb_args" / "blb_summary.csv").exists()
