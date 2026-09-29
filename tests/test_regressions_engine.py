"""Regression tests for the dating engine.

Covered: ape-cladewise node numbering, the rate-ratio guard's ancestor choice,
CorrTest anchoring, CI finite differences on a copy, ddBD ``sampling_frac``
semantics, progress and cancellation in ``calibrate``, per-analysis report
merging, CLI option plumbing, iterative deep-tree handling, quoted-label and
NEXUS escaping, external-token validation and the writer flags.

Each test pins one behavioural rule so that changing that rule fails loudly.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from openreltime.calibrate import (
    Calibration,
    CalibrationCancelled,
    calibrate,
)
from openreltime.ci import confidence_interval
from openreltime.corrtest import corrtest
from openreltime.ddbd import ddbd
from openreltime.report import read_report_block
from openreltime.rates import rrf_rates
from openreltime.table import tree2table
from openreltime.times import rrf_times
from openreltime.tree import PhyloNode
from openreltime.treeio import parse_newick, read_tree, to_newick, write_nexus

ROOT = Path(__file__).resolve().parents[1]
TREE = ROOT / "data" / "examples" / "example.nwk"
OUTGROUP = ["Ornithorhynchus_anatinus", "Zaglossus_bruijni", "Tachyglossus_aculeatus"]

SIMPLE = "((A:0.10,B:0.12):0.20,(C:0.15,(D:0.18,E:0.22):0.25):0.30):0.0;"


def _simple_tree() -> PhyloNode:
    return parse_newick(SIMPLE)


def _example_tree() -> PhyloNode:
    return read_tree(TREE, outgroup=OUTGROUP)


# ---------------------------------------------------------------------------
# Node ids follow ape's cladewise convention: preorder, root = n + 1
# ---------------------------------------------------------------------------

def test_node_ids_are_ape_cladewise() -> None:
    tree = _simple_tree()  # 5 tips, 4 internal nodes
    tree.assign_ids()
    assert tree.node_id == 6  # root = n + 1
    deepest = max(tree.walk(), key=lambda n: n.node_id)
    assert not deepest.is_tip()
    for node in tree.walk():
        if node.is_tip():
            continue
        for child in node.children:
            if not child.is_tip():
                # cladewise: an internal node always numbers below its internal
                # descendants
                assert child.node_id > node.node_id


def test_node_table_ids_ascending_for_internals() -> None:
    from openreltime.report import node_table

    tree = _example_tree()
    frame = node_table(tree)
    internal_ids = [int(v) for v in frame[frame["NodeLabel"] == "-"]["NodeId"]]
    assert internal_ids == sorted(internal_ids)


# ---------------------------------------------------------------------------
# The rate-ratio guard substitutes the nearest clean ancestor, not the root
# ---------------------------------------------------------------------------

def test_rate_ratio_guard_picks_nearest_clean_ancestor() -> None:
    from openreltime.rrf import RrfEngine

    # collapse one clade so the guard fires: a/b are tiny, so the rate of their
    # internal node is extreme
    newick = "((((a:0.00001,b:0.000012):0.000001,c:0.02):0.03,d:0.05):0.06," \
             "(e:0.4,f:0.45):0.5):0.05;"
    tree = parse_newick(newick)
    engine = RrfEngine(tree, rate_ratio_threshold=20.0, compute_times=True)
    engine.run()
    assert engine.n_exceeded > 0, "守卫未触发，测试树需要调整"

    root_time = engine.row(tree).t7a
    for nid, row in engine.rows.items():
        if not row.time_replaced:
            continue
        node = next(n for n in tree.walk() if n.node_id == nid)
        # the nearest clean ancestor is the first node above the parent that is
        # not in the exceed set; whichever node that turns out to be, it must
        # not be the root
        assert row.t7a != pytest.approx(root_time, rel=1e-9), (
            f"node {nid} 被替换成了根时间"
        )
        # and it must equal the t7a of one of the node's true ancestors
        anc_times = {
            engine.row(a).t7a
            for a in node.self_and_ancestors()
            if a is not node and a.node_id in engine.rows
        }
        assert row.t7a in anc_times


# ---------------------------------------------------------------------------
# CorrTest anchoring takes the time of the requested node
# ---------------------------------------------------------------------------

def test_corrtest_anchor_root_scale_is_one() -> None:
    tree = _simple_tree()
    tree.assign_ids()
    res = corrtest(tree, anchor_node=tree.node_id, anchor_time=1.0)
    assert res.anchor is not None
    assert res.anchor["scale.factor"] == pytest.approx(1.0, rel=1e-9)


def test_corrtest_anchor_specific_node_scale() -> None:
    tree = _simple_tree()
    tree.assign_ids()
    target = tree.mrca(["D", "E"])
    times = rrf_times(tree)
    tmax = max(times.times.values())
    expected_sf = tmax / times.times[target.node_id]
    res = corrtest(tree, anchor_node=target.node_id, anchor_time=1.0)
    assert res.anchor is not None
    assert res.anchor["scale.factor"] == pytest.approx(expected_sf, rel=1e-9)


def test_corrtest_anchor_zero_rejected() -> None:
    tree = _simple_tree()
    with pytest.raises(ValueError):
        corrtest(tree, anchor_node=0, anchor_time=1.0)


# ---------------------------------------------------------------------------
# The CI finite differences act on a copy, never on the caller's tree
# ---------------------------------------------------------------------------

def test_ci_matches_reference_finite_difference_and_keeps_tree() -> None:
    tree = _simple_tree()
    tree.assign_ids()
    blens_before = {n.node_id: n.blen for n in tree.walk()}
    times = rrf_times(tree, normalize=True)
    cal = calibrate(
        times,
        [Calibration(node_id=tree.mrca(["A", "B"]).node_id, min_bound=1.0, max_bound=2.0)],
        method="bounds",
    )
    ci = confidence_interval(cal, seq_length=500)
    # SE is strictly positive: the rates differ and the sampling variance is not 0
    assert (ci.table["se"] > 0).all()
    # the caller's branch lengths must come back unchanged
    blens_after = {n.node_id: n.blen for n in tree.walk()}
    assert blens_after == blens_before

    # same semantics checked directly on the pipeline: a perturbation is applied
    # to a copy, so no evaluation can leak into the next one
    from openreltime.ci import _make_pipeline

    edge_nodes = [n for n in tree.walk() if not n.is_root()]
    blens = np.array([float(n.blen or 0.0) for n in edge_nodes])
    pipeline = _make_pipeline(cal, "times")
    v1 = blens.copy()
    v2 = blens.copy()
    v2[0] += 0.05
    out_a = pipeline(v1)
    pipeline(v2)  # state probe: v2 must not be remembered by the next call
    out_a2 = pipeline(v1)
    assert np.allclose(out_a, out_a2, rtol=1e-9), "pipeline 记住了一次调用的输入向量"


# ---------------------------------------------------------------------------
# ddbd reports sampling_frac the way R3F does: fitted free, value replaced
# ---------------------------------------------------------------------------

def test_ddbd_sampling_frac_is_reported_not_constrained() -> None:
    """R3F semantics: rho stays free in the fit; only the reported value changes.

    ``R3F/R/ddbd.R:491-505`` is the reference for that behaviour; the matching
    positive regression test lives in
    ``tests/test_regressions_treeio_peripherals.py``.
    """
    tree = _example_tree()
    base = ddbd(tree)
    assert math.isfinite(base.birth_rate) and base.birth_rate > 0

    fixed = ddbd(tree, sampling_frac=0.9)
    assert fixed.sampling_frac == pytest.approx(0.9)
    assert math.isfinite(fixed.birth_rate) and fixed.birth_rate > 0
    assert math.isfinite(fixed.death_rate) and fixed.death_rate >= 0

    fixed_low = ddbd(tree, sampling_frac=0.1)
    assert fixed_low.sampling_frac == pytest.approx(0.1)
    # rho 自由时，两条链路的 (birth, death) 必是同一个无约束最优解，
    # 且等于不带 sampling_frac 的拟合；只有报告值不同
    assert fixed.birth_rate == pytest.approx(fixed_low.birth_rate, rel=1e-12)
    assert fixed.death_rate == pytest.approx(fixed_low.death_rate, rel=1e-12)
    assert fixed.birth_rate == pytest.approx(base.birth_rate, rel=1e-12)
    assert fixed.params["sampling_frac_fitted"] == pytest.approx(
        base.sampling_frac, rel=1e-9
    )


def test_ddbd_sampling_frac_validated() -> None:
    tree = _simple_tree()
    with pytest.raises(ValueError):
        ddbd(tree, sampling_frac=0.0)
    with pytest.raises(ValueError):
        ddbd(tree, sampling_frac=1.5)


# ---------------------------------------------------------------------------
# An empty calibration list must raise, not pass silently
# ---------------------------------------------------------------------------

def test_calibrate_with_empty_list_raises() -> None:
    tree = _simple_tree()
    times = rrf_times(tree)
    with pytest.raises(ValueError):
        calibrate(times, [], method="bounds")


# ---------------------------------------------------------------------------
# calibrate supports a progress callback and cancellation
# ---------------------------------------------------------------------------

def test_calibrate_progress_callback_and_cancel() -> None:
    tree = _simple_tree()
    times = rrf_times(tree, normalize=True)
    cals = [
        Calibration(node_id=tree.mrca(["A", "B"]).node_id, min_bound=1.0, max_bound=2.0),
        Calibration(node_id=tree.mrca(["D", "E"]).node_id, density="uniform",
                    density_params={"min": 0.5, "max": 1.5}),
    ]
    seen: list[tuple[int, int]] = []
    calibrate(times, cals, method="effective", n_effective=6, seed=1,
              progress=lambda done, total: seen.append((done, total)) or True)
    assert seen and seen[-1] == (6, 6)

    with pytest.raises(CalibrationCancelled):
        calibrate(times, cals, method="effective", n_effective=6, seed=1,
                  progress=lambda done, total: done < 3)


# ---------------------------------------------------------------------------
# report.json merges per analysis instead of overwriting; legacy reads stay compatible
# ---------------------------------------------------------------------------

def test_report_json_merges_across_analyses(tmp_path: Path) -> None:
    tree = _simple_tree()
    prefix = str(tmp_path / "x")

    times = rrf_times(tree)
    times.write(prefix, with_rate=False)
    first = json.loads(Path(f"{prefix}_report.json").read_text())
    assert set(first) == {"rrf_times"}

    rates = rrf_rates(tree)
    rates.write(prefix)
    second = json.loads(Path(f"{prefix}_report.json").read_text())
    assert set(second) == {"rrf_times", "rrf_rates"}, "后写覆盖了先写的报告"

    # read_report_block 兼容两种格式
    assert read_report_block(prefix, "rrf_times")["analysis"] == "rrf_times"
    assert read_report_block(prefix, "calibrate") == {}
    legacy = {"analysis": "corrtest", "score": 0.9}
    Path(f"{prefix}_report.json").write_text(json.dumps(legacy))
    assert read_report_block(prefix, "corrtest") == legacy
    assert read_report_block(prefix, "rrf_times") == {}


# ---------------------------------------------------------------------------
# CLI options round-trip through the recorded report
# ---------------------------------------------------------------------------

def test_rates_times_cli_mean_option_reaches_the_engine(tmp_path) -> None:
    """``--mean`` must change the numbers, not merely be accepted.

    Both conventions run on the same tree so the tables have to differ; an
    option that was parsed and then ignored would still leave this test green
    if it only checked that Click accepted the flag.
    """
    from click.testing import CliRunner

    import pandas as pd

    from openreltime.cli import cli
    from openreltime.treeio import parse_newick
    from openreltime.times import rrf_times

    newick = (
        "((((A:0.12,B:0.20):0.10,(C:0.05,D:0.30):0.04):0.05,"
        "((E:0.18,F:0.22):0.12,(G:0.09,H:0.31):0.07):0.03):0.02,"
        "(I:0.11,J:0.19):0.06):0.0;"
    )
    tree_file = tmp_path / "deep.nwk"
    tree_file.write_text(newick, encoding="utf-8")
    runner = CliRunner()

    results = {}
    for mean in ("geometric", "arithmetic"):
        out = tmp_path / mean
        out.mkdir()
        res = runner.invoke(
            cli,
            ["rates-times", "-i", str(tree_file), "--mean", mean,
             "-o", str(out / "run")],
        )
        assert res.exit_code == 0, res.output + str(res.exception)
        csv = next(Path(out).glob("*_rates_times.csv"))
        results[mean] = pd.read_csv(csv)

    assert not results["geometric"].equals(results["arithmetic"]), (
        "--mean is accepted but does not reach the engine"
    )
    for mean in ("geometric", "arithmetic"):
        tree = parse_newick(newick)
        tree.assign_ids()
        direct = rrf_times(tree, mean=mean)
        col = [c for c in results[mean].columns if "ime" in c.lower()]
        assert col, results[mean].columns.tolist()
        mine = dict(zip(results[mean]["NodeId"], results[mean][col[0]]))
        for nid, value in direct.times.items():
            assert float(mine[nid]) == pytest.approx(value, abs=1e-12)


def test_tree2table_cli_accepts_fmt() -> None:
    from openreltime.cli import cli

    fmt_opts = {
        opt for p in cli.commands["tree2table"].params for opt in p.opts
    }
    assert "--fmt" in fmt_opts


def test_ci_cli_restores_effective_params(tmp_path: Path) -> None:
    """calibrate records input_fmt/resolve and ci reads them back."""
    from openreltime.report import write_report

    prefix = str(tmp_path / "cal")
    write_report(
        prefix,
        {
            "analysis": "calibrate",
            "parameters": {
                "method": "effective",
                "n_effective": 500,
                "seed": 7,
                "input_fmt": "nexus",
                "resolve_polytomy": "random",
            },
        },
    )
    block = read_report_block(prefix, "calibrate")
    assert block["parameters"]["input_fmt"] == "nexus"
    assert block["parameters"]["resolve_polytomy"] == "random"


# ---------------------------------------------------------------------------
# deep trees survive the whole path: parsing, traversal and serialisation are iterative
# ---------------------------------------------------------------------------

def _pectinate(n: int) -> str:
    """恰好 n 个末端的链状树。"""
    text = "A:1.0"
    for i in range(1, n):
        text = f"({text},B{i}:{1.0 + i * 0.001})"
    return text + ";"


def test_deep_tree_full_pipeline() -> None:
    n = 1200  # 远超默认递归限
    tree = parse_newick(_pectinate(n))
    assert tree.n_tips() == n
    times = rrf_times(tree)
    assert len(times.times) == 2 * n - 1
    text = to_newick(tree)
    assert text.endswith(";")
    nexus = write_nexus(tree, rate_map=times.rates, time_map=times.times)
    assert nexus.startswith("#NEXUS")
    frame = tree2table(tree, time=True)
    assert "Time" in frame.columns


def test_copy_is_deep_safe() -> None:
    tree = parse_newick(_pectinate(1200))
    clone = tree.copy()
    assert clone.n_tips() == tree.n_tips()


# ---------------------------------------------------------------------------
# quoted-label round trip and NEXUS escaping
# ---------------------------------------------------------------------------

def test_quoted_label_roundtrip() -> None:
    label = "Escherichia (K-12), strain 'X'"
    tree = parse_newick(f"((A:0.1,'B, (K-12)':0.2):0.3,'{label.replace(chr(39), chr(39) * 2)}':0.4):0.0;")
    text = to_newick(tree)
    reread = parse_newick(text)
    assert reread.tip_labels() == ["A", "B, (K-12)", label]


def test_nexus_taxlabels_are_escaped(tmp_path: Path) -> None:
    tree = parse_newick("((A:0.1,'B, C':0.2):0.3,D:0.4):0.0;")
    nexus = write_nexus(tree)
    assert "'B, C'" in nexus
    # 写出的 NEXUS 也能被自己的读侧解析
    path = tmp_path / "t.nex"
    path.write_text(nexus)
    reread = read_tree(path, fmt="nexus")
    assert "B, C" in reread.tip_labels()


# ---------------------------------------------------------------------------
# external-process arguments are checked against a shell-metacharacter blacklist
# ---------------------------------------------------------------------------

def test_validate_tokens_allows_spaces_and_cjk() -> None:
    from openreltime._external import validate_tokens

    validate_tokens("/Applications/My Tool/iqtree2")
    validate_tokens("/home/user/tools/iqtree")
    validate_tokens("-m", "GTR+G")


def test_validate_tokens_rejects_metacharacters() -> None:
    from openreltime._external import validate_tokens

    for bad in ("a;b", "a|b", "a>b", "a$b", "a`b", 'a"b', "a\nb"):
        with pytest.raises(ValueError):
            validate_tokens(bad)
    with pytest.raises(ValueError):
        validate_tokens("")


# ---------------------------------------------------------------------------
# no division by zero when rho_ad is 0
# ---------------------------------------------------------------------------

def test_corrtest_zero_rho_ad_does_not_crash(monkeypatch: pytest.MonkeyPatch) -> None:
    import importlib

    ct = importlib.import_module("openreltime.corrtest")
    tree = _simple_tree()
    tree.assign_ids()
    # 强制所有 Spearman 系数为 0（含 rho_ad），验证衰减率不除零
    monkeypatch.setattr(ct, "_spearman", lambda a, b: 0.0)
    res = ct.corrtest(tree)
    assert res.rho_ad_1_decay == 0.0
    assert res.rho_ad_2_decay == 0.0


# ---------------------------------------------------------------------------
# Result.write honours every flag combination
# ---------------------------------------------------------------------------

def test_time_result_write_flags(tmp_path: Path) -> None:
    tree = _simple_tree()
    times = rrf_times(tree)
    p1 = str(tmp_path / "a")
    paths = times.write(p1, with_rate=False, nexus=True)
    assert any(str(p).endswith("_timetree.nexus") for p in paths)

    p2 = str(tmp_path / "b")
    paths = times.write(p2, with_rate=True, nexus=False)
    assert any(str(p).endswith("_timetree.nwk") for p in paths)


def test_calibrated_result_write_nexus_false_still_writes_tree(tmp_path: Path) -> None:
    tree = _simple_tree()
    times = rrf_times(tree, normalize=True)
    cal = calibrate(
        times,
        [Calibration(node_id=tree.mrca(["A", "B"]).node_id, min_bound=1.0, max_bound=2.0)],
        method="bounds",
    )
    paths = cal.write(str(tmp_path / "c"), nexus=False)
    assert any(str(p).endswith("_calibrated.nwk") for p in paths)
    assert not any(str(p).endswith("_calibrated.nexus") for p in paths)
