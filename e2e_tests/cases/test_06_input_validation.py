"""End-to-end tests for input validation and clean error handling."""

from __future__ import annotations

import re
from pathlib import Path


from e2e_tests.conftest import (  # noqa: E402
    OUTGROUP_STR,
    TREE,
)


def _modify_nwk(tmp_out: Path, pattern: str, repl: str) -> Path:
    """Return a modified Newick file under tmp_out."""
    text = TREE.read_text()
    new_text = re.sub(pattern, repl, text, count=1)
    out = tmp_out / "bad.nwk"
    out.write_text(new_text)
    return out


def test_negative_branch_length_rejected(run_cli_raw, tmp_out):
    """Trees with negative branch lengths are rejected."""
    text = TREE.read_text()
    # Modify an ingroup branch length so it survives outgroup pruning
    text = text.replace("Rhyncholestes_raphanurus:0.05537700", "Rhyncholestes_raphanurus:-0.05537700", 1)
    bad = tmp_out / "neg.nwk"
    bad.write_text(text)
    proc = run_cli_raw("times", "-i", str(bad), "--outgroup", OUTGROUP_STR, "-o", "neg")
    assert proc.returncode != 0
    assert "negative" in proc.stderr.lower()


def test_duplicate_tip_rejected(run_cli_raw, tmp_out):
    """Trees with duplicated tip labels are rejected."""
    bad = _modify_nwk(tmp_out, r"Zaglossus_bruijni", "Zaglossus_bruijni,Zaglossus_bruijni")
    proc = run_cli_raw("times", "-i", str(bad), "--outgroup", OUTGROUP_STR, "-o", "dup")
    assert proc.returncode != 0
    assert "duplicate" in proc.stderr.lower()


def test_polytomy_rejected_by_default(run_cli_raw, tmp_out):
    """Polytomies raise by default and can be resolved with --resolve random."""
    text = TREE.read_text()
    # Turn the first ingroup clade into a trichotomy by duplicating one tip with a new name
    text = text.replace(
        "(Rhyncholestes_raphanurus:0.05537700,Caenolestes_fuliginosus:0.05342500):0.05176100",
        "(Rhyncholestes_raphanurus:0.05537700,Caenolestes_fuliginosus:0.05342500,Caenolestes_fuliginosus_copy:0.010000):0.05176100",
    )
    bad = tmp_out / "poly.nwk"
    bad.write_text(text)
    proc = run_cli_raw("times", "-i", str(bad), "--outgroup", OUTGROUP_STR, "-o", "poly_err")
    assert proc.returncode != 0
    assert "polytom" in proc.stderr.lower()

    # With --resolve random it should succeed
    proc2 = run_cli_raw(
        "times", "-i", str(bad), "--outgroup", OUTGROUP_STR,
        "--resolve", "random", "-o", "poly_ok",
    )
    assert proc2.returncode == 0


def test_too_few_ingroup_tips_rejected(run_cli_raw, tmp_out):
    """An outgroup that is not a complete clade is rejected."""
    proc = run_cli_raw(
        "times",
        "-i", str(TREE),
        "--outgroup",
        "Ornithorhynchus_anatinus,Zaglossus_bruijni,Tachyglossus_aculeatus,Rhyncholestes_raphanurus",
        "-o", "few",
    )
    assert proc.returncode != 0
    assert "complete clade" in proc.stderr.lower()


def test_incomplete_outgroup_rejected(run_cli_raw, tmp_out):
    """A multi-tip outgroup that is not a complete clade errors by default."""
    proc = run_cli_raw(
        "times",
        "-i", str(TREE),
        "--outgroup", "Ornithorhynchus_anatinus,Zaglossus_bruijni",
        "-o", "inc_og",
    )
    assert proc.returncode != 0
    assert "complete clade" in proc.stderr.lower()

    # With warn it should succeed
    proc2 = run_cli_raw(
        "times",
        "-i", str(TREE),
        "--outgroup", "Ornithorhynchus_anatinus,Zaglossus_bruijni",
        "--outgroup-check", "warn",
        "-o", "inc_og_warn",
    )
    assert proc2.returncode == 0


def test_empty_calibration_row_rejected(run_cli_raw, tmp_out):
    """A calibration row with neither bound nor density is rejected."""
    bad = tmp_out / "empty_cal.tsv"
    bad.write_text(
        "node_id\ttaxon_set\tmin_bound\tmax_bound\tdensity\tdensity_params\n"
        "-\tHomo_sapiens|Pan_troglodytes\t.\t.\t.\t.\n"
    )
    proc = run_cli_raw(
        "calibrate",
        "-i", str(TREE),
        "-c", str(bad),
        "--outgroup", OUTGROUP_STR,
        "-o", "empty_cal",
    )
    assert proc.returncode != 0
    assert "bound" in proc.stderr.lower()


def test_density_only_with_bounds_method_rejected(run_cli_raw, tmp_out):
    """Density-only calibration with --method bounds is rejected."""
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
        "-o", "dens_bounds",
    )
    assert proc.returncode != 0
    assert "effective" in proc.stderr.lower()


def test_missing_input_file(run_cli_raw):
    """A missing tree file gives a clean CLI error."""
    proc = run_cli_raw("times", "-i", "does_not_exist.nwk", "-o", "missing")
    assert proc.returncode != 0
    assert "not found" in proc.stderr.lower() or "no such" in proc.stderr.lower()
