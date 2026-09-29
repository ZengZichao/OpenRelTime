"""Core RRF engine: local postorder estimation + global preorder adjustment.

Implements the relative rate framework (RRF) analytical solution of

    Tamura K, Tao Q, Kumar S.  *The relative rate framework reveals the
    timings of speciation events.*  Mol Biol Evol 35:1770-1782 (2018),
    doi:10.1093/molbev/msy044

following the behaviour of the GPL-3 reference implementation R3F
(Tao & Kumar 2025, github.com/articledoten/R3F).  Where the paper leaves
implementation details open, R3F semantics were replicated exactly and are
cited inline (``R3F/R/rrf_times.R`` line numbers; see also
``THIRD-PARTY-NOTICES.md``).

Two averaging conventions are available:

* ``"geometric"`` (default) -- the R3F/msy044 eqs 28-42 convention;
* ``"arithmetic"`` -- the original PNAS 2012 semantics (msy044 eqs 1-27),
  provided for cross-validation; it is *not* bit-compatible with R3F.

Both conventions apply to the *whole* recursion: the local solve **and** the
folding of child-subtree lengths into the parent's composite lengths (``l1 =
sqrt(l1*l2) + l5`` geometric vs. ``l1 = (l1+l2)/2 + l5`` arithmetic; msy044
fig. 2 legend: "La = b5 + 1/2(b1 + b2) ... when using the arithmetic mean").
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Optional

from . import _constants as C
from .tree import PhyloNode

logger = logging.getLogger("openreltime")

__all__ = ["RrfRow", "RrfEngine"]


@dataclass
class RrfRow:
    """Per-node RRF workspace, mirroring the ``RRF.mat`` matrix of R3F.

    ``r5``/``r6`` are the local rates of the two child lineages of the node
    (relative to the node's own local root); ``t7`` is the local root age of
    the node's clade; ``r5a``/``r6a``/``t7a`` are the globally adjusted
    quantities.  ``l1..l6`` are the composite lengths used in the local
    solution (``l1 = fold(l1a, l1b) + l5`` folding, cf. msy044 figs. 2-3;
    ``fold`` is :func:`RrfEngine._fold`, geometric in R3F and arithmetic under
    the PNAS 2012 convention).
    """

    node_id: int
    child1: Optional[PhyloNode] = None
    child2: Optional[PhyloNode] = None
    l1: float = 0.0
    l2: float = 0.0
    l3: float = 0.0
    l4: float = 0.0
    l5: float = 0.0
    l6: float = 0.0
    r5: Optional[float] = None
    r6: Optional[float] = None
    t7: Optional[float] = None
    r5a: Optional[float] = None
    r6a: Optional[float] = None
    t7a: Optional[float] = None
    #: set when the rate-ratio guard replaced the node time
    time_replaced: bool = False
    #: set when a tip-grandchild rate was backfilled into r5/r6
    backfilled: list[str] = field(default_factory=list)


class RrfEngine:
    """Postorder + preorder RRF computation for one rooted binary tree.

    Parameters
    ----------
    tree:
        Rooted binary ingroup tree (see :func:`openreltime.read_tree`).
    mean:
        ``"geometric"`` (R3F-compatible) or ``"arithmetic"`` (PNAS 2012);
        the choice applies to the local solve *and* to the composite-length
        folding (:meth:`_fold`).
    rate_ratio_threshold:
        ``None`` (the R3F ``rrf_rates`` default,
        :data:`openreltime._constants.RATE_RATIO_GUARD_DISABLED`) disables the
        guard; a number enables it (R3F ``rrf_times`` default,
        :data:`openreltime._constants.RATE_RATIO_THRESHOLD`).  The guard only
        ever rewrites node *times* (``rrf_times.R`` lines 353-373), so it is
        inert when ``compute_times`` is false.
    compute_times:
        Compute local/global node times (``rrf_times`` pipeline).
    """

    def __init__(
        self,
        tree: PhyloNode,
        *,
        mean: str = "geometric",
        rate_ratio_threshold: Optional[float] = C.RATE_RATIO_GUARD_DISABLED,
        compute_times: bool = True,
    ) -> None:
        if mean not in ("geometric", "arithmetic"):
            raise ValueError(f"mean must be 'geometric' or 'arithmetic', got {mean!r}")
        # ids are the row keys of the workspace; ensure they exist even when
        # the tree was built with parse_newick directly
        tree.assign_ids()
        self.tree = tree
        self.mean = mean
        self.rate_ratio_threshold = rate_ratio_threshold
        self.compute_times = compute_times
        self.rows: dict[int, RrfRow] = {}
        self.n_exceeded: int = 0
        for node in tree.walk():
            if not node.is_tip():
                self.rows[node.node_id] = RrfRow(node_id=node.node_id)

    # -- public API ---------------------------------------------------------

    def run(self) -> "RrfEngine":
        self._postorder_local()
        self._preorder_adjust()
        # The guard belongs to the times pass only: R3F rrf_times.R lines
        # 353-373 rewrite ``t7.adjust``, and rrf_rates.R carries no guard at all
        # (rates are never replaced), so a threshold without ``compute_times``
        # has nothing to act on.
        if self.compute_times and self.rate_ratio_threshold is not None:
            self._rate_ratio_guard()
        return self

    def row(self, node: PhyloNode) -> RrfRow:
        return self.rows[node.node_id]

    # -- postorder: local estimation -----------------------------------------

    def _postorder_local(self) -> None:
        """Per-node local solution (R3F rrf_times.R lines 79-273)."""
        tree = self.tree
        for node in tree.iter_postorder():
            if node.is_tip():
                continue
            row = self.row(node)
            c1, c2 = node.children
            row.child1, row.child2 = c1, c2
            if c1.is_tip() and c2.is_tip():
                # 2-clade: record stem lengths only; rates arrive from parent
                row.l5 = c1.blen or 0.0
                row.l6 = c2.blen or 0.0
                continue
            if not c1.is_tip() and not c2.is_tip():
                self._four_clade(node, row, c1, c2)
            elif not c1.is_tip() and c2.is_tip():
                self._three_clade_internal_first(node, row, c1, c2)
            else:
                self._three_clade_tip_first(node, row, c1, c2)

    # -- composite lengths --------------------------------------------------

    def _fold(self, a: float, b: float) -> float:
        """Average two sister lengths of one clade under the chosen ``mean``.

        Geometric reproduces R3F exactly (``rrf_times.R`` lines 105-108,
        178-179, 230-231); the arithmetic branch is the PNAS 2012 / msy044
        eqs 1-27 convention, whose composite depths are "averaged branch
        lengths", ``La = b5 + 1/2(b1 + b2)`` (msy044 fig. 2 legend).  Using
        the geometric fold here while solving locally with the arithmetic
        mean mixes the two conventions and breaks ``b = r * t`` on trees
        deeper than two levels.
        """
        if self.mean == "geometric":
            return math.sqrt(a * b)
        return (a + b) / 2

    # composite lengths of an internal child: l_a = fold(l1, l2) + l5 etc.
    def _composite(self, row: RrfRow) -> tuple[float, float]:
        return self._fold(row.l1, row.l2) + row.l5, self._fold(row.l3, row.l4) + row.l6

    def _four_clade(
        self, node: PhyloNode, row: RrfRow, c1: PhyloNode, c2: PhyloNode
    ) -> None:
        r1row, r2row = self.row(c1), self.row(c2)
        row.l1, row.l2 = self._composite(r1row)
        row.l3, row.l4 = self._composite(r2row)
        row.l5 = c1.blen or 0.0
        row.l6 = c2.blen or 0.0
        l1, l2, l3, l4, l5, l6 = row.l1, row.l2, row.l3, row.l4, row.l5, row.l6
        # epsilon guards, exactly as placed in R3F rrf_times.R lines 114-131
        if l1 == 0:
            l1 = C.EPS_R3F
        if l2 == 0:
            l2 = C.EPS_R3F
        if l3 == 0:
            l3 = C.EPS_R3F
        if l4 == 0:
            l4 = C.EPS_R3F
        if self._fold(l3, l4) + l6 == 0:
            l6 = C.EPS_R3F_SMALL
        if self._fold(l1, l2) + l5 == 0:
            l5 = C.EPS_R3F_SMALL
        la = self._fold(l1, l2) + l5
        lb = self._fold(l3, l4) + l6
        if self.mean == "geometric":
            r5 = math.sqrt(la) / math.sqrt(lb)
            r6 = math.sqrt(lb) / math.sqrt(la)
            row.t7 = math.sqrt(la) * math.sqrt(lb)
        else:
            # msy044 eqs 19-27, written on the arithmetic composite lengths so
            # that the local identity la = r5 * t7 holds by construction
            r5 = 2 * la / (la + lb)
            r6 = 2 * lb / (la + lb)
            row.t7 = (la + lb) / 2
        row.r5, row.r6 = r5, r6
        # grandchild rates (used only when a grandchild is a tip)
        if self.mean == "geometric":
            rr1 = math.sqrt(l1) * math.sqrt(la) / (math.sqrt(l2) * math.sqrt(lb))
            rr2 = math.sqrt(l2) * math.sqrt(la) / (math.sqrt(l1) * math.sqrt(lb))
            rr3 = math.sqrt(l3) * math.sqrt(lb) / (math.sqrt(l4) * math.sqrt(la))
            rr4 = math.sqrt(l4) * math.sqrt(lb) / (math.sqrt(l3) * math.sqrt(la))
            t_back1 = math.sqrt(l1 * l2) * math.sqrt(lb) / math.sqrt(la)
            t_back2 = math.sqrt(l3 * l4) * math.sqrt(la) / math.sqrt(lb)
        else:
            rr1 = 4 * l1 * la / ((l1 + l2) * (la + lb))
            rr2 = 4 * l2 * la / ((l1 + l2) * (la + lb))
            rr3 = 4 * l3 * lb / ((l3 + l4) * (la + lb))
            rr4 = 4 * l4 * lb / ((l3 + l4) * (la + lb))
            t_back1 = (l1 + l2) * (la + lb) / (4 * la)
            t_back2 = (l3 + l4) * (la + lb) / (4 * lb)
        if c1.children[0].is_tip():
            r1row.r5 = rr1
            if self.compute_times:
                r1row.t7 = t_back1
            r1row.backfilled.append("r5")
        if c1.children[1].is_tip():
            r1row.r6 = rr2
            if self.compute_times:
                r1row.t7 = t_back1
            r1row.backfilled.append("r6")
        if c2.children[0].is_tip():
            r2row.r5 = rr3
            if self.compute_times:
                r2row.t7 = t_back2
            r2row.backfilled.append("r5")
        if c2.children[1].is_tip():
            r2row.r6 = rr4
            if self.compute_times:
                r2row.t7 = t_back2
            r2row.backfilled.append("r6")

    def _three_clade_internal_first(
        self, node: PhyloNode, row: RrfRow, c1: PhyloNode, c2: PhyloNode
    ) -> None:
        """3-clade with ``c1`` internal and ``c2`` a tip."""
        r1row = self.row(c1)
        row.l1, row.l2 = self._composite(r1row)
        row.l3 = 0.0
        row.l4 = 0.0
        row.l5 = c1.blen or 0.0
        row.l6 = c2.blen or 0.0
        l1, l2, l5, l6 = row.l1, row.l2, row.l5, row.l6
        # epsilon guards as in R3F rrf_times.R lines 187-195
        if l1 == 0:
            l1 = C.EPS_R3F
        if l2 == 0:
            l2 = C.EPS_R3F
        if l6 == 0:
            l6 = C.EPS_R3F
        la = self._fold(l1, l2) + l5
        if self.mean == "geometric":
            r5 = math.sqrt(la) / math.sqrt(l6)
            r6 = math.sqrt(l6) / math.sqrt(la)
            row.t7 = math.sqrt(la) * math.sqrt(l6)
            rr1 = math.sqrt(l1) * math.sqrt(la) / math.sqrt(l2 * l6)
            rr2 = math.sqrt(l2) * math.sqrt(la) / math.sqrt(l1 * l6)
            t_back = math.sqrt(l1 * l2) * math.sqrt(l6) / math.sqrt(la)
        else:
            # msy044 eqs 1-18 on the arithmetic composite lengths (the tip
            # lineage is its own composite: lb = fold(0, 0) + l6 = l6)
            r5 = 2 * la / (la + l6)
            r6 = 2 * l6 / (la + l6)
            row.t7 = (la + l6) / 2
            rr1 = 4 * l1 * la / ((l1 + l2) * (la + l6))
            rr2 = 4 * l2 * la / ((l1 + l2) * (la + l6))
            t_back = (l1 + l2) * (la + l6) / (4 * la)
        row.r5, row.r6 = r5, r6
        if c1.children[0].is_tip():
            r1row.r5 = rr1
            if self.compute_times:
                r1row.t7 = t_back
            r1row.backfilled.append("r5")
        if c1.children[1].is_tip():
            r1row.r6 = rr2
            if self.compute_times:
                r1row.t7 = t_back
            r1row.backfilled.append("r6")

    def _three_clade_tip_first(
        self, node: PhyloNode, row: RrfRow, c1: PhyloNode, c2: PhyloNode
    ) -> None:
        """3-clade with ``c1`` a tip and ``c2`` internal."""
        r2row = self.row(c2)
        row.l1 = 0.0
        row.l2 = 0.0
        row.l3, row.l4 = self._composite(r2row)
        row.l5 = c1.blen or 0.0
        row.l6 = c2.blen or 0.0
        l3, l4, l5, l6 = row.l3, row.l4, row.l5, row.l6
        # epsilon guards as in R3F rrf_times.R lines 237-245
        if l5 == 0:
            l5 = C.EPS_R3F
        if l3 == 0:
            l3 = C.EPS_R3F
        if l4 == 0:
            l4 = C.EPS_R3F
        lb = self._fold(l3, l4) + l6
        if self.mean == "geometric":
            r5 = math.sqrt(l5) / math.sqrt(lb)
            r6 = math.sqrt(lb) / math.sqrt(l5)
            row.t7 = math.sqrt(l5) * math.sqrt(lb)
            rr3 = math.sqrt(l3) * math.sqrt(lb) / math.sqrt(l4 * l5)
            rr4 = math.sqrt(l4) * math.sqrt(lb) / math.sqrt(l3 * l5)
            t_back = math.sqrt(l3 * l4) * math.sqrt(l5) / math.sqrt(lb)
        else:
            # msy044 eqs 1-18 (the tip lineage is its own composite:
            # la = fold(0, 0) + l5 = l5)
            r5 = 2 * l5 / (l5 + lb)
            r6 = 2 * lb / (l5 + lb)
            row.t7 = (l5 + lb) / 2
            rr3 = 4 * l3 * lb / ((l3 + l4) * (l5 + lb))
            rr4 = 4 * l4 * lb / ((l3 + l4) * (l5 + lb))
            t_back = (l3 + l4) * (l5 + lb) / (4 * lb)
        row.r5, row.r6 = r5, r6
        if c2.children[0].is_tip():
            r2row.r5 = rr3
            if self.compute_times:
                r2row.t7 = t_back
            r2row.backfilled.append("r5")
        if c2.children[1].is_tip():
            r2row.r6 = rr4
            if self.compute_times:
                r2row.t7 = t_back
            r2row.backfilled.append("r6")

    # -- preorder: global adjustment ------------------------------------------

    def _preorder_adjust(self) -> None:
        """Multiply local rates by ancestral rates (rrf_times.R lines 282-351).

        The "grandparent" special case is replicated verbatim: when a node
        has a tip child, the corresponding backfilled slot already refers to
        the parent's local frame, so the *grandparent's* adjusted rate is
        used instead.  The asymmetry of the R3F reference (the no-grandparent
        case refreshes ``t7a`` only for the first-child slot) is kept as is.
        """
        tree = self.tree
        root_row = self.row(tree)
        root_row.r5a = root_row.r5
        root_row.r6a = root_row.r6
        if self.compute_times:
            root_row.t7a = root_row.t7
        for node in tree.iter_preorder():
            if node.is_root() or node.is_tip():
                continue
            parent = node.parent
            assert parent is not None
            prow = self.row(parent)
            nrow = self.row(node)
            is_first = parent.children[0] is node
            r_anc = prow.r5a if is_first else prow.r6a
            if r_anc is None:
                r_anc = 1.0
            nrow.r5a = None if nrow.r5 is None else nrow.r5 * r_anc
            nrow.r6a = None if nrow.r6 is None else nrow.r6 * r_anc
            if self.compute_times and nrow.t7 is not None:
                nrow.t7a = nrow.t7 / r_anc
            grandpa = parent.parent
            for slot, child_idx in ((5, 0), (6, 1)):
                if not node.children[child_idx].is_tip():
                    continue
                if grandpa is None:
                    # parent is the root (R3F lines 298-300 / 309-311 and
                    # 327-329 / 338-340): only the r5 branch refreshes t7a
                    if slot == 5 and self.compute_times and nrow.t7 is not None:
                        nrow.t7a = nrow.t7
                else:
                    gprow = self.row(grandpa)
                    gp_first = grandpa.children[0] is parent
                    r_g = gprow.r5a if gp_first else gprow.r6a
                    if r_g is None:
                        r_g = 1.0
                    if slot == 5:
                        nrow.r5a = None if nrow.r5 is None else nrow.r5 * r_g
                    else:
                        nrow.r6a = None if nrow.r6 is None else nrow.r6 * r_g
                    if self.compute_times and nrow.t7 is not None:
                        nrow.t7a = nrow.t7 / r_g

    # -- rate-ratio guard -------------------------------------------------------

    def _rate_ratio_guard(self) -> None:
        """Replace extreme node times with the nearest clean ancestor's time.

        Replicates R3F rrf_times.R lines 353-373, where the reference hard-codes
        ``rate.ratio = 20`` (:data:`openreltime._constants.RATE_RATIO_THRESHOLD`):
        the ancestor with the largest ``ape`` id outside the exceeding set
        wins — under cladewise (pre-order) numbering that is the ancestor
        *closest* to the replaced node.  Only ``t7.adjust`` is rewritten;
        rates are untouched, which is why the rates-only pass has no guard.
        """
        threshold = float(self.rate_ratio_threshold)
        exceed_r5: list[int] = []
        exceed_r6: list[int] = []
        for nid, row in self.rows.items():
            if row.r5a is not None and (row.r5a > threshold or row.r5a < 1 / threshold):
                exceed_r5.append(nid)
            if row.r6a is not None and (row.r6a > threshold or row.r6a < 1 / threshold):
                exceed_r6.append(nid)
        exceed_node_ids: set[int] = set()
        for nid in exceed_r5:
            child = self.rows[nid].child1
            if child is not None and not child.is_tip():
                exceed_node_ids.add(child.node_id)
        for nid in exceed_r6:
            child = self.rows[nid].child2
            if child is not None and not child.is_tip():
                exceed_node_ids.add(child.node_id)

        for nid in exceed_r5:
            self._replace_time(self.rows[nid].child1, exceed_node_ids)
        for nid in exceed_r6:
            self._replace_time(self.rows[nid].child2, exceed_node_ids)
        self.n_exceeded = len(exceed_node_ids)
        if exceed_node_ids:
            logger.warning(
                "%d node time(s) replaced by an ancestral time because the "
                "adjusted rate ratio exceeded %.0f",
                len(exceed_node_ids),
                threshold,
            )

    def _replace_time(self, child: Optional[PhyloNode], exceed_node_ids: set[int]) -> None:
        if child is None or child.is_tip():
            return
        ancestors = [
            anc
            for anc in child.self_and_ancestors()
            if anc is not child and anc.node_id in self.rows
        ]
        clean = [a for a in ancestors if a.node_id not in exceed_node_ids]
        if not clean:
            return
        # R3F takes max(ape id): under cladewise numbering the ancestors of a
        # node have descending ids toward the root, so max() is the ancestor
        # closest to *child* (the R3F "nearest clean ancestor" semantics)
        chosen = max(clean, key=lambda a: a.node_id)
        crow = self.rows[chosen.node_id]
        self.rows[child.node_id].t7a = crow.t7a
        self.rows[child.node_id].time_replaced = True
