"""Monophyly assessment and outgroup complete-clade checking (G8 gate).

Covers the two outgroup input forms (single tip / multiple tips), the
multi-tip complete-clade (MRCA completeness) requirement, group-name
outgroups and calibrations resolved from a user taxon table or
auto-detected tip-label prefixes, and the non-monophyly warnings.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest
from click.testing import CliRunner

from openreltime import Calibration, calibrate
from openreltime.cli import cli
from openreltime.taxonomy import (
    check_monophyly,
    detect_groups_from_labels,
    parse_taxon_table,
)
from openreltime.times import rrf_times
from openreltime.treeio import parse_newick, read_tree, root_outgroup

# Mus/Rattus/Homo/Pan/Elephas/Loxodonta, 8 tips, rooted at 0.0
TREE = (
    "(((Mus_musculus:0.1,Mus_spretus:0.2):0.3,"
    "(Rattus_norvegicus:0.4,Rattus_rattus:0.5):0.6):0.7,"
    "((Homo_sapiens:0.2,Pan_troglodytes:0.3):0.4,"
    "(Elephas_maximus:0.5,Loxodonta_africana:0.6):0.7):0.8):0.0;"
)


@pytest.fixture()
def tree():
    t = parse_newick(TREE)
    t.assign_ids()
    return t


# ---------------------------------------------------------------------------
# check_monophyly
# ---------------------------------------------------------------------------

def test_monophyletic_group_detected(tree) -> None:
    res = check_monophyly(tree, ["Mus_musculus", "Mus_spretus"], name="Mus")
    assert res.is_monophyletic
    assert res.extra_tips == ()
    assert res.missing == ()


def test_non_monophyletic_group_lists_extra_tips(tree) -> None:
    # elephant pair + human: their MRCA is the (Homo,Pan)/(Elephas,Lox) clade,
    # so chimpanzee also descends from it
    res = check_monophyly(
        tree, ["Elephas_maximus", "Loxodonta_africana", "Homo_sapiens"]
    )
    assert not res.is_monophyletic
    assert set(res.extra_tips) == {"Pan_troglodytes"}


def test_unknown_names_reported_as_missing(tree) -> None:
    res = check_monophyly(tree, ["Mus_musculus", "Nofoundus_nowhere"])
    assert not res.is_monophyletic
    assert res.missing == ("Nofoundus_nowhere",)


def test_single_tip_group_is_monophyletic(tree) -> None:
    assert check_monophyly(tree, ["Homo_sapiens"]).is_monophyletic


# ---------------------------------------------------------------------------
# outgroup forms (requirement 1)
# ---------------------------------------------------------------------------

def test_outgroup_single_tip_and_multi_tip_complete_clade(tmp_path: Path) -> None:
    path = tmp_path / "t.nwk"
    path.write_text(TREE)
    one = read_tree(path, outgroup=["Homo_sapiens"])
    assert set(one.tip_labels()) == {
        "Mus_musculus",
        "Mus_spretus",
        "Rattus_norvegicus",
        "Rattus_rattus",
        "Pan_troglodytes",
        "Elephas_maximus",
        "Loxodonta_africana",
    }
    two = read_tree(path, outgroup=["Mus_musculus", "Mus_spretus"])
    assert "Mus_musculus" not in two.tip_labels()
    assert two.is_binary()


def test_incomplete_multi_tip_outgroup_rejected(tmp_path: Path) -> None:
    path = tmp_path / "t.nwk"
    path.write_text(TREE)
    with pytest.raises(ValueError, match="not a complete clade"):
        read_tree(path, outgroup=["Mus_musculus", "Rattus_norvegicus"])


def test_incomplete_outgroup_warn_mode_proceeds(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    path = tmp_path / "t.nwk"
    path.write_text(TREE)
    with caplog.at_level(logging.WARNING, logger="openreltime"):
        t = read_tree(
            path,
            outgroup=["Mus_musculus", "Rattus_norvegicus"],
            outgroup_check="warn",
        )
    assert any("not a complete clade" in r.message for r in caplog.records)
    # both named outgroup tips are still removed
    assert "Mus_musculus" not in t.tip_labels()
    assert "Rattus_norvegicus" not in t.tip_labels()


def test_root_outgroup_check_utility(tree) -> None:
    pruned = root_outgroup(tree, ["Mus_musculus"])
    assert "Mus_musculus" not in pruned.tip_labels()
    assert len(pruned.tips()) == len(tree.tips()) - 1
    from openreltime.treeio import check_outgroup_clade

    assert check_outgroup_clade(tree, ["Mus_musculus"]) == []
    assert check_outgroup_clade(
        tree, ["Mus_musculus", "Mus_spretus", "Rattus_norvegicus"]
    ) == ["Rattus_rattus"]


# ---------------------------------------------------------------------------
# group-name outgroups via auto-detection and taxon table (requirement 2)
# ---------------------------------------------------------------------------

def test_outgroup_by_group_name_auto_detected(tmp_path: Path) -> None:
    path = tmp_path / "t.nwk"
    path.write_text(TREE)
    # "Mus" is not a tip; it is auto-detected as the genus of Mus_* tips
    t = read_tree(path, outgroup=["Mus"])
    assert set(t.tip_labels()) == {
        "Rattus_norvegicus",
        "Rattus_rattus",
        "Homo_sapiens",
        "Pan_troglodytes",
        "Elephas_maximus",
        "Loxodonta_africana",
    }


def test_outgroup_by_group_name_from_taxon_table(tmp_path: Path) -> None:
    path = tmp_path / "t.nwk"
    table_path = tmp_path / "taxa.tsv"
    path.write_text(TREE)
    table_path.write_text(
        "tip\tgroup\n"
        "Mus_musculus\tMurinae\n"
        "Mus_spretus\tMurinae\n"
        "Rattus_norvegicus\tMurinae\n"
        "Rattus_rattus\tMurinae\n"
        "Homo_sapiens\tPrimates\n"
        "Pan_troglodytes\tPrimates\n"
        "Elephas_maximus\tAfrotheria\n"
        "Loxodonta_africana\tAfrotheria\n"
    )
    t = read_tree(path, outgroup=["Murinae"], taxon_table=parse_taxon_table(table_path))
    assert set(t.tip_labels()) == {
        "Homo_sapiens",
        "Pan_troglodytes",
        "Elephas_maximus",
        "Loxodonta_africana",
    }


def test_outgroup_by_non_monophyletic_group_name_warns_then_errors(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A group judged non-monophyletic raises a warning (requirement 2) and,
    by default, the outgroup completeness check then rejects it (req. 1)."""
    path = tmp_path / "t.nwk"
    table_path = tmp_path / "taxa.tsv"
    path.write_text(TREE)
    table_path.write_text(
        "tip\tgroup\n"
        "Mus_musculus\tBroken\n"
        "Mus_spretus\tBroken\n"
        "Rattus_norvegicus\tBroken\n"  # Rattus_rattus deliberately absent
        "Homo_sapiens\tOthers\n"
        "Pan_troglodytes\tOthers\n"
        "Elephas_maximus\tOthers\n"
        "Loxodonta_africana\tOthers\n"
    )
    with caplog.at_level(logging.WARNING, logger="openreltime"):
        with pytest.raises(ValueError, match="not a complete clade"):
            read_tree(
                path,
                outgroup=["Broken"],
                taxon_table=parse_taxon_table(table_path),
            )
    assert any("NOT monophyletic" in r.message for r in caplog.records)


def test_mixed_tip_and_group_outgroup_tokens(tmp_path: Path) -> None:
    path = tmp_path / "t.nwk"
    path.write_text(TREE)
    # "Mus" expands to the two Mus_* tips; together with both Rattus tips the
    # outgroup covers the complete rodent clade
    t = read_tree(path, outgroup=["Mus", "Rattus_norvegicus", "Rattus_rattus"])
    assert set(t.tip_labels()) == {
        "Homo_sapiens",
        "Pan_troglodytes",
        "Elephas_maximus",
        "Loxodonta_africana",
    }


# ---------------------------------------------------------------------------
# taxonomy table and label auto-detection
# ---------------------------------------------------------------------------

def test_parse_taxon_table_aliases_and_csv(tmp_path: Path) -> None:
    csv_path = tmp_path / "taxa.csv"
    csv_path.write_text(
        "label,classification\n"
        "Mus_musculus,Murinae\n"
        "Homo_sapiens,Primates\n"
    )
    table = parse_taxon_table(csv_path)
    assert table == {"Murinae": ["Mus_musculus"], "Primates": ["Homo_sapiens"]}

    # two columns with unrecognised headers: (first -> tip, second -> group)
    fallback = tmp_path / "plain.tsv"
    fallback.write_text("sample\tclade_x\nMus_musculus\tMurinae\nHomo_sapiens\tPrimates\n")
    table2 = parse_taxon_table(fallback)
    assert table2 == {"Murinae": ["Mus_musculus"], "Primates": ["Homo_sapiens"]}


def test_parse_taxon_table_duplicate_tip_rejected(tmp_path: Path) -> None:
    path = tmp_path / "dup.tsv"
    path.write_text("tip\tgroup\nA\tX\nA\tY\n")
    with pytest.raises(ValueError, match="twice"):
        parse_taxon_table(path)


def test_detect_groups_from_labels() -> None:
    groups = detect_groups_from_labels(
        ["Homo_sapiens", "Homo_heidelbergensis", "Pan_paniscus|ABC1", "Singleton"]
    )
    assert groups["Homo"] == ["Homo_sapiens", "Homo_heidelbergensis"]
    assert groups["Pan"] == ["Pan_paniscus|ABC1"]
    assert groups["Singleton"] == ["Singleton"]


# ---------------------------------------------------------------------------
# calibrations by group name (requirement 2)
# ---------------------------------------------------------------------------

def test_calibration_by_group_name(tmp_path: Path) -> None:
    path = tmp_path / "t.nwk"
    path.write_text(TREE)
    t = read_tree(path)
    times = rrf_times(t)
    cal = calibrate(
        times,
        [Calibration(taxon_set=frozenset(["Afrotheria"]), min_bound=0.5, max_bound=2.0)],
        taxon_table={
            "Afrotheria": ["Elephas_maximus", "Loxodonta_africana"],
            "Others": ["Mus_musculus", "Rattus_norvegicus"],
        },
    )
    record = cal.calibrations[0]
    assert set(record.target) == {"Elephas_maximus", "Loxodonta_africana"}


def test_calibration_non_monophyletic_group_warns(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    path = tmp_path / "t.nwk"
    path.write_text(TREE)
    t = read_tree(path)
    times = rrf_times(t)
    with caplog.at_level(logging.WARNING, logger="openreltime"):
        cal = calibrate(
            times,
            [Calibration(taxon_set=frozenset(["Mus_musculus", "Homo_sapiens"]), min_bound=0.01)],
        )
    assert any("NOT monophyletic" in r.message for r in caplog.records)
    # calibration still attaches to the MRCA of the requested tips
    assert len(cal.calibrations) == 1


# ---------------------------------------------------------------------------
# CLI monophyly subcommand
# ---------------------------------------------------------------------------

def test_cli_monophyly_reports_verdicts(tmp_path: Path) -> None:
    path = tmp_path / "t.nwk"
    path.write_text(TREE)
    runner = CliRunner()
    # auto-detected groups (genus prefixes)
    result = runner.invoke(
        cli,
        ["monophyly", "-i", str(path), "--groups", "Mus,Homo", "-o", str(tmp_path / "chk")],
    )
    assert result.exit_code == 0, result.output
    import pandas as pd

    frame = pd.read_csv(tmp_path / "chk_monophyly.csv")
    verdicts = dict(zip(frame["group"], frame["is_monophyletic"]))
    assert bool(verdicts["Mus"]) is True
    assert bool(verdicts["Homo"]) is True

    # table-based group that is not monophyletic
    table_path = tmp_path / "taxa.tsv"
    table_path.write_text(
        "tip\tgroup\n"
        "Mus_musculus\tOdd\n"
        "Mus_spretus\tOdd\n"
        "Homo_sapiens\tOdd\n"  # MRCA of Odd is the whole tree
        "Rattus_norvegicus\tRest\n"
        "Rattus_rattus\tRest\n"
        "Pan_troglodytes\tRest\n"
        "Elephas_maximus\tRest\n"
        "Loxodonta_africana\tRest\n"
    )
    result2 = runner.invoke(
        cli,
        [
            "monophyly",
            "-i",
            str(path),
            "--taxon-table",
            str(table_path),
            "--groups",
            "Odd",
            "-o",
            str(tmp_path / "chk2"),
        ],
    )
    assert result2.exit_code == 0, result2.output
    frame2 = pd.read_csv(tmp_path / "chk2_monophyly.csv")
    assert not bool(frame2["is_monophyletic"].iloc[0])
    extra = frame2["extra_tips"].iloc[0]
    assert "Pan_troglodytes" in str(extra) and "Rattus_rattus" in str(extra)


def test_cli_monophyly_warns_on_stderr(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """The monophyly subcommand still succeeds but raises a warning for the
    non-monophyletic set (logging goes to the stderr stream)."""
    path = tmp_path / "t.nwk"
    path.write_text(TREE)
    runner = CliRunner()
    with caplog.at_level(logging.WARNING, logger="openreltime"):
        result = runner.invoke(
            cli,
            [
                "monophyly",
                "-i",
                str(path),
                "--tips",
                "Mus_musculus,Homo_sapiens",
                "-o",
                str(tmp_path / "chk"),
            ],
        )
    assert result.exit_code == 0, result.output
    assert any("NOT monophyletic" in r.message for r in caplog.records)
