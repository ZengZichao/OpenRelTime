"""G2a/G2b gates: regression against stored R3F golden outputs.

The golden files in ``data/golden/r3f/`` were generated with R3F
(GPL-3, R 4.5.3, ape 5.8.1) on the bundled 274-tip mammalian tree
(``data/examples/example.nwk``, dos Reis et al. 2012); the generation
script and versions are recorded in ``data/PROVENANCE.md``.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from compare_golden import compare  # noqa: E402

from openreltime.times import rrf_times  # noqa: E402
from openreltime.treeio import read_tree, to_newick  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TREE = ROOT / "data" / "examples" / "example.nwk"
GOLDEN = ROOT / "data" / "golden" / "r3f"
OUTGROUP = ["Ornithorhynchus_anatinus", "Zaglossus_bruijni", "Tachyglossus_aculeatus"]


@pytest.fixture(scope="module")
def comparison() -> dict[str, float]:
    return compare(TREE, GOLDEN, OUTGROUP, with_rate=True)


@pytest.fixture(scope="module")
def tree():
    return read_tree(TREE, outgroup=OUTGROUP)


def test_g2a_rate_regression(comparison: dict[str, float]) -> None:
    """G2a: slope >= 0.999 and median relative deviation <= 1e-6."""
    assert comparison["rate_slope"] >= 0.999
    assert comparison["rate_median_rel"] <= 1e-6


def test_g2b_time_regression(comparison: dict[str, float]) -> None:
    """G2b: slope >= 0.999 and all clade times matched."""
    assert comparison["n_common"] == 541  # 271 tips + 270 internal nodes
    assert comparison["time_slope"] >= 0.999
    assert comparison["time_max_rel"] <= 1e-6


def test_rooted_topology_matches_ape(tree) -> None:
    """Rooting + outgroup pruning reproduces ape::root + drop.tip exactly."""
    expected = (GOLDEN / "example_rooted_dropped.nwk").read_text().strip()
    mine = to_newick(tree)
    assert _sorted_pairs(mine) == _sorted_pairs(expected)


def _sorted_pairs(newick: str) -> list[tuple[str, float]]:
    """Canonicalise a tree as sorted (tip label, root-to-tip path) pairs."""
    from openreltime.treeio import parse_newick

    tree = parse_newick(newick if newick.endswith(";") else newick + ";")

    def depths(n, acc: float = 0.0):
        out: list[tuple[str, float]] = []
        if n.is_tip():
            out.append((n.label or "", round(acc, 9)))
        for c in n.children:
            out.extend(depths(c, acc + (c.blen or 0.0)))
        return out

    return sorted(depths(tree))


def test_rates_times_table_matches_golden(tree) -> None:
    """The rates-times table equals R3F's _RRF_table.csv node by node."""
    res = rrf_times(tree)
    golden = pd.read_csv(GOLDEN / "example_RRF_table.csv")

    def by_clade(df: pd.DataFrame, value: str) -> dict:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from compare_golden import clade_map_from_csv

        if "Des1" in df.columns and df["NodeLabel"].astype(str).eq("-").any():
            clades = (
                clade_map_from_csv(df)
                if "Des1" in df.columns
                else {}
            )
        else:
            clades = {int(r.NodeId): frozenset([r.NodeLabel]) for r in df.itertuples()}
        out = {}
        for r in df.itertuples():
            nid = int(r.NodeId)
            if nid in clades:
                out[clades[nid]] = float(getattr(r, value))
        return out

    g_time = by_clade(golden, "Time")
    m_time = {n.clade: res.times[n.node_id] for n in tree.walk()}
    common = set(g_time) & set(m_time)
    diffs = [abs(g_time[c] - m_time[c]) for c in common if g_time[c] > 0]
    assert max(diffs) < 1e-6


def test_tree2table_matches_golden(tmp_path: Path) -> None:
    """tree2table reproduces R3F's branch-length pairs, matched by clade."""
    golden = pd.read_csv(GOLDEN / "example_tree2table.csv", dtype=str)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from compare_golden import clade_map_from_csv

    from openreltime.treeio import read_tree as _read

    raw_tree = _read(TREE)  # golden tree2table ran on the unpruned tree
    gclades = clade_map_from_csv(golden)
    golden_pairs: dict[frozenset, tuple[float, float]] = {}
    for r in golden.itertuples():
        nid = int(r.NodeId)
        if r.NodeLabel != "-" or str(nid) not in {str(k) for k in gclades}:
            continue
        clade = gclades[nid]
        golden_pairs[clade] = (float(r.Brlen1), float(r.Brlen2))

    mine_pairs: dict[frozenset, tuple[float, float]] = {}
    for n in raw_tree.walk():
        if n.is_tip():
            continue
        mine_pairs[n.clade] = (
            float(n.children[0].blen or 0.0),
            float(n.children[1].blen or 0.0),
        )
    assert set(golden_pairs) == set(mine_pairs)
    for clade, (g1, g2) in golden_pairs.items():
        m1, m2 = mine_pairs[clade]
        pair_golden = {round(g1, 6), round(g2, 6)}
        pair_mine = {round(m1, 6), round(m2, 6)}
        assert pair_golden == pair_mine, f"clade mismatch: {sorted(clade)[:2]}"


def test_corrttest_matches_r3f_score(tree) -> None:
    """CorrTest on the golden tree reproduces the stored R3F CorrScore.

    The golden value ``score = 0.9996`` is recorded in
    ``data/golden/r3f/example_sr0_corrtest.txt``.  A P-value band plus
    ``score > 0.9`` would not detect an internally shifted rate vector, because
    mis-paired internal-node rates still leave the tip-anchored Spearman
    correlations above 0.9, so the score itself is pinned to the reference value
    alongside the band.
    """
    from openreltime.corrtest import corrtest
    from openreltime.rates import rrf_rates

    golden_line = (GOLDEN / "example_sr0_corrtest.txt").read_text().splitlines()
    golden_score = float(
        next(ln for ln in golden_line if ln.startswith("score")).split("=")[1]
    )
    res = corrtest(tree)
    assert res.p_band == "P-value < 0.001"
    assert abs(res.score - golden_score) < 1e-4  # the printed value is 5 digits

    # The printed golden value carries only five significant digits, so it cannot
    # serve as an error budget.  Compare against the unrounded double instead:
    # OpenRelTime reproduces it exactly.
    full = (GOLDEN / "example_corrtest_full_precision.txt").read_text()
    exact = float([l for l in full.splitlines()
                   if l.startswith("sr0")][0].split("\t")[1])
    assert res.score == exact, (res.score, exact)

    # the per-node rates CorrTest consumes must be the same physical quantity
    # the engine reports, not a vector shifted by one generation
    from openreltime.corrtest import _rrf_rate_vectors

    ids, rates, _d1, _d2 = _rrf_rate_vectors(tree)
    engine_rates = rrf_rates(tree, rate_ratio_threshold=None).rates
    internal = [n for n in tree.walk() if not n.is_tip()]
    assert internal  # guard against a vacuous comparison
    for node in internal:
        assert engine_rates[node.node_id] == pytest.approx(
            rates[ids.index(node.node_id)], abs=1e-9
        ), f"internal node {node.node_id} rate differs"


def test_corrttest_lag_decay_features_are_not_degenerate(tree) -> None:
    """lag-2/lag-3 decay features must be computed from real 2-/3-step pairs.

    Pairing on the (child, grandchild) generation gap alone drops the start
    node, which makes ``rho_ad_1_decay`` an exact structural zero rather than a
    measured decay; the assertion therefore distinguishes "not computed" from
    "computed and small".
    """
    from openreltime.corrtest import corrtest

    res = corrtest(tree)
    assert res.rho_ad != 0.0
    assert math.isfinite(res.rho_ad_1_decay) and math.isfinite(res.rho_ad_2_decay)
    # a genuine multi-generation decay must not be an exact structural zero
    assert not (res.rho_ad_1_decay == 0.0 and res.rho_ad_2_decay == 0.0)


def test_ddbd_close_to_r3f(tree) -> None:
    """ddBD parameters match the R3F run within optimizer tolerance."""
    from openreltime.ddbd import ddbd

    # anchor_node is what switches anchoring on in R3F (ddbd.R:393-396); with no
    # anchor_node the fit stays in relative time, which is what this file records.
    res = ddbd(tree)
    golden_line = (GOLDEN / "example_noanchor_ddbd.txt").read_text().splitlines()[1]
    gb, gd, gs = (float(v) for v in golden_line.split("\t"))
    assert res.birth_rate == pytest.approx(gb, rel=5e-4)
    assert res.death_rate == pytest.approx(gd, rel=5e-4)
    assert res.sampling_frac == pytest.approx(gs, rel=5e-3)


def test_remaining_golden_files_are_reproduced(tree) -> None:
    """Every stored R3F golden file is consumed by at least one assertion.

    Otherwise the anchored ddBD path, the resampled CorrTest value, the timed
    trees and the rates-only comparison would ship as reference data that no
    test reads.
    """
    from openreltime.corrtest import corrtest
    from openreltime.ddbd import ddbd

    # anchored ddBD: R3F's own worked example (anchor.node = 272, anchor.time = 1.85)
    anchor_golden = [
        float(v) for v in (GOLDEN / "example_anchor_ddbd.txt").read_text()
        .splitlines()[1].split("\t")
    ]
    ad = ddbd(tree, anchor_node=272, anchor_time=1.85)
    assert ad.birth_rate == pytest.approx(anchor_golden[0], rel=5e-4)
    assert ad.death_rate == pytest.approx(anchor_golden[1], rel=5e-4)
    assert ad.sampling_frac == pytest.approx(anchor_golden[2], rel=5e-3)
    assert ad.scale_factor == pytest.approx(1.85, rel=1e-9)

    # resampled CorrTest golden value
    sr50 = [
        corrtest(tree, sister_resample=50, seed=s).score for s in (1, 2, 3)
    ]
    for score in sr50:
        assert abs(score - 0.99959) < 1e-4, score

    # timed trees: Newick heights and NEXUS rate annotations
    gold_times = pd.read_csv(GOLDEN / "example_RRF_times.csv")
    root_time = float(gold_times["Time"].max())
    tt = read_tree(GOLDEN / "example_RRF_timetree.nwk")
    assert len(tt.tips()) == 271
    depths = []
    for tip in tt.tips():
        node, acc = tip, 0.0
        while node is not None:
            acc += float(node.blen or 0.0)
            node = node.parent
        depths.append(acc)
    assert max(depths) == pytest.approx(root_time, abs=1e-6)

    nexus_text = (GOLDEN / "example_RRF_timetree.nexus").read_text()
    assert "[&rate=" in nexus_text, "golden NEXUS lost its rate annotations"
    # the golden NEXUS opens with a hol-style "[#NEXUS ...]" comment block, which
    # the suffix-free auto-detector does not claim; pass the format explicitly.
    back = read_tree(GOLDEN / "example_RRF_timetree.nexus", fmt="nexus")
    assert len(back.tips()) == 271

    rates_nwk = (GOLDEN / "example_RRF_rates.nwk").read_text()
    assert rates_nwk.count(":") > 271  # every lineage carries a rate-length

    # the rates-only branch of the comparison helper is exercised too
    no_rate = compare(TREE, GOLDEN, OUTGROUP, with_rate=False)
    assert "rate_slope" not in no_rate
    assert no_rate["n_common"] == 541
