"""Convert a tree to an R3F-compatible node table (``tree2table``)."""

from __future__ import annotations

import pandas as pd

from .report import node_table
from .tree import PhyloNode

__all__ = ["tree2table"]


def tree2table(tree: PhyloNode, *, time: bool = False) -> pd.DataFrame:
    """Convert a tree to a table (R3F ``tree2table`` equivalent).

    Parameters
    ----------
    tree:
        Parsed tree (see :func:`openreltime.read_tree`).
    time:
        ``False`` (default) reports branch lengths; ``True`` reports for
        each internal node the longest cumulative branch length from the
        node down to its tips in a ``Time`` column (identical to the ape
        ``branching.times`` height on ultrametric trees; on
        non-ultrametric trees the two definitions differ — this column is
        the max-depth-to-tips one).
    """
    df = node_table(tree)
    brlen1: list[str] = ["-"] * len(df)
    brlen2: list[str] = ["-"] * len(df)
    child1 = {n.node_id: n.children[0] for n in tree.walk() if not n.is_tip()}
    child2 = {n.node_id: n.children[1] for n in tree.walk() if not n.is_tip()}
    for i, row in df.iterrows():
        if row["NodeLabel"] == "-":
            nid = int(row["NodeId"])
            brlen1[i] = repr(child1[nid].blen)
            brlen2[i] = repr(child2[nid].blen)
    df["Brlen1"] = brlen1
    df["Brlen2"] = brlen2
    if time:
        # node ages: longest distance from the node to any of its tips,
        # computed iteratively in post-order (deep-tree safe)
        depth: dict[int, float] = {}
        for n in tree.iter_postorder():
            if n.is_tip():
                continue
            depth[n.node_id] = max(
                (c.blen or 0.0) + depth.get(c.node_id, 0.0) for c in n.children
            )
        time_col = ["-"] * len(df)
        for i, row in df.iterrows():
            if row["NodeLabel"] == "-":
                time_col[i] = repr(depth[int(row["NodeId"])])
        df["Time"] = time_col
    return df
