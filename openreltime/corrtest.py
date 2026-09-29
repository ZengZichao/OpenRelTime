"""CorrTest: detection of autocorrelated evolutionary rates.

Implements the fixed-coefficient logistic model of

    Tao Q, Tamura K, Mello B, Kumar S.  *Bioinformatic identification of
    autocorrelated rates in molecular evolution.*  Mol Biol Evol 36:811-824
    (2019), doi:10.1093/molbev/msz014

following the GPL-3 R3F implementation (``R3F/R/corrtest.R``); the
normalisation constants and logistic coefficients were extracted from lines
555-583 and are stored in :mod:`openreltime._constants`.
"""

from __future__ import annotations

import logging
import math
from typing import Any, Optional

import numpy as np
from scipy.stats import spearmanr

from . import _constants as C
from .report import CorrTestResult
from .rrf import RrfEngine
from .tree import PhyloNode

logger = logging.getLogger("openreltime")

__all__ = ["corrtest"]


def _rrf_rate_vectors(
    tree: PhyloNode, threshold: Optional[float] = None
):
    """Local RRF rates per node exactly as prepared in corrtest.R lines 359-368.

    Returns ``(node_ids, rates, des1, des2)`` where ``node_ids`` are 1..n+nnode
    in ape order (tips first), ``rates`` are the adjusted rates with root = 1
    and ``des1``/``des2`` are 0 for tips.

    The rate attached to a node is always the *adjusted rate of the lineage
    leading to it*, i.e. it is stored in the parent row's ``r5a``/``r6a`` slot.
    ``corrtest.R:360-361`` makes this explicit::

        nodeID    <- c(RRF.mat[,2], RRF.mat[,3], tips.num+1)
        RRF.rates <- c(RRF.mat[,"r5.adjust"], RRF.mat[,"r6.adjust"], 1)

    so tips and internal nodes are treated identically, the only special case
    being the root, whose rate is 1 by construction.
    """
    if threshold is None:
        # Source: openreltime/_constants.RATE_RATIO_THRESHOLD (rrf_times.R:354)
        threshold = C.RATE_RATIO_THRESHOLD
    engine = RrfEngine(tree, mean="geometric", compute_times=True)
    engine.run()
    # R3F zeroes out local rates above the threshold *before* adjustment
    for row in engine.rows.values():
        if row.r5 is not None and row.r5 > threshold:
            row.r5 = 0.0
        if row.r6 is not None and row.r6 > threshold:
            row.r6 = 0.0
    # redo the adjustment with zeroed locals (corrtest.R lines 289-353)
    root_row = engine.row(tree)
    root_row.r5a = root_row.r5
    root_row.r6a = root_row.r6
    for node in tree.iter_preorder():
        if node.is_root() or node.is_tip():
            continue
        parent = node.parent
        assert parent is not None
        prow = engine.row(parent)
        nrow = engine.row(node)
        r_anc = prow.r5a if parent.children[0] is node else prow.r6a
        if r_anc is None:
            r_anc = 1.0
        nrow.r5a = None if nrow.r5 is None else nrow.r5 * r_anc
        nrow.r6a = None if nrow.r6 is None else nrow.r6 * r_anc
        for slot, child_idx in ((5, 0), (6, 1)):
            if not node.children[child_idx].is_tip():
                continue
            grandpa = parent.parent
            if grandpa is None:
                continue
            gprow = engine.row(grandpa)
            r_g = gprow.r5a if grandpa.children[0] is parent else gprow.r6a
            if r_g is None:
                r_g = 1.0
            if slot == 5:
                nrow.r5a = None if nrow.r5 is None else nrow.r5 * r_g
            else:
                nrow.r6a = None if nrow.r6 is None else nrow.r6 * r_g

    # ape/cladewise id order: tips 1..n (in tip order), internal n+1.. in
    # pre-order (root = n + 1); iterate ids ascending to reproduce it
    id_to_node = {n.node_id: n for n in tree.walk()}
    node_ids: list[int] = []
    rates: list[float] = []
    des1: list[int] = []
    des2: list[int] = []
    for nid in sorted(id_to_node):
        node = id_to_node[nid]
        node_ids.append(nid)
        if node.is_tip():
            # tip rate: stored in the parent's adjusted child slot
            parent = node.parent
            prow = engine.row(parent)
            value = prow.r5a if parent.children[0] is node else prow.r6a
            rates.append(float(value) if value is not None else 0.0)
            des1.append(0)
            des2.append(0)
        else:
            # internal node: its own rate is the adjusted lineage rate stored
            # in the *parent's* slot, exactly as for tips.  ``row(node).r5a``
            # would instead give the rate of the node's first child, i.e. the
            # whole vector shifted down one generation.
            if node.is_root():
                # root rate is 1 by definition; its adjusted slots carry the
                # rates of its two children, which are picked up as their own
                # parent-slot values below
                rates.append(1.0)
            else:
                parent = node.parent
                prow = engine.row(parent)
                value = prow.r5a if parent.children[0] is node else prow.r6a
                rates.append(float(value) if value is not None else 0.0)
            des1.append(node.children[0].node_id)
            des2.append(node.children[1].node_id)
    return node_ids, rates, des1, des2


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    rho = spearmanr(a, b).statistic
    if rho is None or math.isnan(rho):
        return 0.0
    return float(rho)


def corrtest(
    tree: PhyloNode,
    *,
    sister_resample: int = 0,
    seed: Optional[int] = None,
    anchor_node: Optional[int] = None,
    anchor_time: float = 0.0,
    **params: Any,
) -> CorrTestResult:
    """Test the hypothesis of independence of evolutionary rates.

    Parameters
    ----------
    tree:
        Rooted binary ingroup tree with branch lengths.
    sister_resample:
        Number of random sister-pair swaps (use > 50 for trees with fewer
        than ~50 tips); ``0`` disables resampling.
    seed:
        Seed for the sister-pair resampling.
    anchor_node, anchor_time:
        When given (a node id and its absolute time), relative rates are
        converted to absolute rates and their mean/SD reported.

    Returns
    -------
    CorrTestResult
        CorrScore in [0, 1] and its P-value band.
    """
    node_ids, rates_arr, des1, des2 = _rrf_rate_vectors(tree)
    rates = np.array(rates_arr, dtype=float)
    rates[~np.isfinite(rates)] = 0.0
    id_index = {nid: i for i, nid in enumerate(node_ids)}

    internal = [i for i in range(len(node_ids)) if des1[i] != 0]

    # ---- sister-lineage correlation ---------------------------------------
    pairs1 = np.array([rates[id_index[des1[i]]] for i in internal])
    pairs2 = np.array([rates[id_index[des2[i]]] for i in internal])
    keep = (pairs1 != 0) & (pairs2 != 0)
    r1, r2 = pairs1[keep], pairs2[keep]
    if len(r1) < 3:
        raise ValueError(
            "tree has fewer than 3 rate-bearing sister pairs; CorrTest needs "
            "a larger tree"
        )
    rho_s = _spearman(r1, r2)

    if sister_resample:
        rng = np.random.default_rng(seed)
        rhos = []
        for _ in range(int(sister_resample)):
            swap = rng.integers(0, 2, size=len(r1))
            r1s = np.concatenate([r1[swap == 0], r2[swap == 1]])
            r2s = np.concatenate([r2[swap == 0], r1[swap == 1]])
            rhos.append(_spearman(r1s, r2s))
        rho_s = float(np.mean(rhos))

    # ---- ancestor-descendant correlation -----------------------------------
    # Everything below is built in *node-id space* so that the two entries per
    # internal node, the parent lookup and the lag descent follow corrtest.R
    # lines 425-540 literally.  Node-id space is also what keeps the lag
    # features apart: positions in ``internal``, positions in ``node_ids`` and
    # node ids are three incompatible index spaces, and conflating them advances
    # every start node by the same number of generations.  Lag-2 and lag-3 then
    # degenerate into near-copies of the lag-1 statistic instead of pairing a
    # node with its grandchildren / great-grandchildren while holding the
    # *original* ancestor's rate fixed.
    def _parent_id(nid: int) -> int:
        for j in internal:
            if nid == des1[j] or nid == des2[j]:
                return node_ids[j]
        return 0

    nodes: list[int] = []
    anc_nodes: list[int] = []
    anc_rates: list[float] = []
    des_nodes: list[int] = []
    des_rates: list[float] = []
    for i in internal:
        anc = _parent_id(node_ids[i])
        for slot in (des1[i], des2[i]):
            nodes.append(node_ids[i])
            anc_nodes.append(anc)
            anc_rates.append(rates[i])
            des_nodes.append(slot)
            des_rates.append(rates[id_index[slot]])

    def _filtered_correlation(
        a: list[float], b: list[float], has_parent: list[bool]
    ) -> float:
        # corrtest.R: drop entries whose start node has no parent (the root),
        # then drop the zero-rate placeholders left by the guard
        ra = [x for x, y, k in zip(a, b, has_parent) if k and x != 0 and y != 0]
        rd = [y for x, y, k in zip(a, b, has_parent) if k and x != 0 and y != 0]
        if not ra:
            return 0.0
        return _spearman(np.array(ra), np.array(rd))

    has_parent = [k != 0 for k in anc_nodes]
    rho_ad = _filtered_correlation(anc_rates, des_rates, has_parent)
    rho_ad_all = [rho_ad]

    for lag in (2, 3):
        lag_anc: list[float] = []
        lag_des: list[float] = []
        lag_has_parent: list[bool] = []
        for i in range(len(nodes)):
            # the two entries whose start node is the selected descendant;
            # empty when that descendant is a tip -> corrtest.R's `next`
            group = [q for q, nv in enumerate(nodes) if nv == des_nodes[i]]
            if not group:
                continue
            for _ in range(lag - 2):
                descended: list[int] = []
                for q in group:
                    descended.extend(
                        r for r, nv in enumerate(nodes) if nv == des_nodes[q]
                    )
                group = descended
                if not group:
                    break
            if not group:
                continue
            for q in group:
                # the ancestor is always the original start node, replicated
                # once per descendant (corrtest.R: rep(anc.rates[i], ...))
                lag_anc.append(anc_rates[i])
                lag_des.append(des_rates[q])
                lag_has_parent.append(anc_nodes[i] != 0)
        if not lag_anc:
            rho_ad_all.append(0.0)
            continue
        rho_ad_all.append(
            _filtered_correlation(lag_anc, lag_des, lag_has_parent)
        )

    if rho_ad_all[0] == 0:
        logger.warning(
            "ancestor-descendant correlation is 0 (degenerate rate data); "
            "the lag-2/lag-3 decays are set to 0"
        )
        rho_ad_1_decay = 0.0
        rho_ad_2_decay = 0.0
    else:
        rho_ad_1_decay = (rho_ad_all[1] - rho_ad_all[0]) / rho_ad_all[0]
        rho_ad_2_decay = (rho_ad_all[2] - rho_ad_all[0]) / rho_ad_all[0]

    # ---- fixed-coefficient logistic model (corrtest.R lines 555-567) -------
    rho_s_norm = (rho_s - C.CORRTEST_NORM_MEAN[0]) / C.CORRTEST_NORM_SD[0]
    rho_ad_norm = (rho_ad_all[0] - C.CORRTEST_NORM_MEAN[1]) / C.CORRTEST_NORM_SD[1]
    d1_norm = (rho_ad_1_decay - C.CORRTEST_NORM_MEAN[2]) / C.CORRTEST_NORM_SD[2]
    d2_norm = (rho_ad_2_decay - C.CORRTEST_NORM_MEAN[3]) / C.CORRTEST_NORM_SD[3]
    b0, b1, b2, b3, b4 = C.CORRTEST_COEFFS
    score = 1.0 / (
        1.0 + math.exp(-(b0 + b1 * rho_s_norm + b2 * rho_ad_norm + b3 * d1_norm + b4 * d2_norm))
    )
    p_band = C.corrtest_p_band(score)

    anchor_stats: Optional[dict[str, float]] = None
    if anchor_node is not None:
        if anchor_node <= 0:
            raise ValueError(
                f"anchor node id must be a positive internal-node id, got "
                f"{anchor_node}"
            )
        engine = RrfEngine(tree, compute_times=True)
        engine.run()
        tvals = np.array(
            [engine.row(n).t7a or 0.0 for n in tree.walk() if not n.is_tip()]
        )
        rel_time = tvals / tvals.max()
        rel_time[rel_time < 0] = 0.0
        # anchor node lookup: walk() is pre-order, so index by the node's
        # *position in the same array*, never by id arithmetic
        id_to_pos = {
            n.node_id: i
            for i, n in enumerate(
                n for n in tree.walk() if not n.is_tip()
            )
        }
        target = next(
            (n for n in tree.walk() if not n.is_tip() and n.node_id == anchor_node),
            None,
        )
        if target is None:
            raise KeyError(f"anchor node id {anchor_node} not found among internal nodes")
        if rel_time[id_to_pos[anchor_node]] <= 0:
            raise ValueError(
                f"anchor node {anchor_node} has zero relative time and cannot "
                f"be used for rescaling"
            )
        sf = anchor_time / rel_time[id_to_pos[anchor_node]]
        rel_rates = np.array(
            [
                engine.row(n).r5a or 0.0
                for n in tree.walk()
                if not n.is_tip()
            ]
            + [
                engine.row(n).r6a or 0.0
                for n in tree.walk()
                if not n.is_tip()
            ]
        )
        abs_rates = rel_rates / sf
        anchor_stats = {
            "mean.rate": float(np.mean(abs_rates)),
            "sd.rate": float(np.std(abs_rates, ddof=1)),
            "scale.factor": float(sf),
        }

    merged = {
        "sister_resample": sister_resample,
        "seed": seed,
        **params,
    }
    return CorrTestResult(
        rho_s=float(rho_s),
        rho_ad=float(rho_ad_all[0]),
        rho_ad_1_decay=float(rho_ad_1_decay),
        rho_ad_2_decay=float(rho_ad_2_decay),
        score=float(score),
        p_band=p_band,
        anchor=anchor_stats,
        params=merged,
    )
