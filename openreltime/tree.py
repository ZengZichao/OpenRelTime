"""Lightweight binary phylogenetic tree used throughout OpenRelTime.

The RRF analytical solution of Tamura et al. (2018, MBE 35:1770-1782) is a
node-local algorithm that needs controllable post-/pre-order traversals,
composite-length caches and clade identities.  Neither Bio.Phylo nor DendroPy
exposed these efficiently (ADR-001/ADR-002, indexed in docs/adr/README.md), hence
this minimal self-contained structure.
"""

from __future__ import annotations

from typing import Iterable, Iterator, Optional


class PhyloNode:
    """A single node of a rooted phylogenetic tree.

    Attributes
    ----------
    label:
        Tip name for leaves, optional internal label otherwise.
    blen:
        Branch length (substitutions/site) of the edge connecting this node
        to its parent.  ``None`` for the root.
    support:
        Optional support value parsed from internal node labels.
    rate:
        Optional relative rate read back from a ``[&rate=...]`` Newick
        annotation (see :func:`openreltime.treeio.read_tree`).
    children:
        Ordered list of child nodes (0, 1 or 2 elements).  The order follows
        the Newick reading order and mirrors ``ape``/``phangorn`` semantics.
    parent:
        Parent node or ``None`` for the root.
    node_id:
        Integer identifier assigned by :meth:`assign_ids`.
    internal:
        ``True`` for a node that represents clade structure, i.e. that was
        read or created as an internal node.  It is set by :meth:`add_child`
        and by the Newick reader, and it is what distinguishes a real tip
        from an internal node that lost all of its children (a "ghost tip",
        removed by :func:`prune_childless_internal_nodes`).
    """

    __slots__ = (
        "label",
        "blen",
        "support",
        "rate",
        "children",
        "parent",
        "node_id",
        "internal",
        "_clade",
        "_tips",
    )

    def __init__(
        self,
        label: Optional[str] = None,
        blen: Optional[float] = None,
        support: Optional[float] = None,
        rate: Optional[float] = None,
        internal: bool = False,
    ) -> None:
        self.label = label
        self.blen = blen
        self.support = support
        self.rate = rate
        self.children: list["PhyloNode"] = []
        self.parent: Optional["PhyloNode"] = None
        self.node_id: int = 0
        self.internal = internal
        self._clade: Optional[frozenset[str]] = None
        self._tips: Optional[list[str]] = None

    # -- structure ---------------------------------------------------------

    def add_child(self, child: "PhyloNode", blen: Optional[float] = None) -> "PhyloNode":
        if blen is not None:
            child.blen = blen
        child.parent = self
        self.children.append(child)
        self.internal = True
        self._invalidate()
        return child

    def _invalidate(self) -> None:
        for node in self.self_and_ancestors():
            node._clade = None
            node._tips = None

    def is_tip(self) -> bool:
        return not self.children

    def is_root(self) -> bool:
        return self.parent is None

    def is_binary(self) -> bool:
        return all(len(n.children) in (0, 2) for n in self.walk())

    # -- traversals --------------------------------------------------------
    # All traversals are iterative: phylogenetic trees can be far deeper than
    # Python's default recursion limit (~1000), and a pectinate tree of a few
    # thousand tips must not raise RecursionError (see docs/methods.md).

    def walk(self) -> Iterator["PhyloNode"]:
        """Yield every node of the subtree in pre-order (root first)."""
        stack = [self]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(reversed(node.children))

    def iter_postorder(self) -> Iterator["PhyloNode"]:
        """Yield nodes children-before-parents (deepest first)."""
        stack = [self]
        order: list["PhyloNode"] = []
        while stack:
            node = stack.pop()
            order.append(node)
            stack.extend(node.children)
        yield from reversed(order)

    def iter_preorder(self) -> Iterator["PhyloNode"]:
        """Yield nodes parents-before-children (root first)."""
        return self.walk()

    def self_and_ancestors(self) -> Iterator["PhyloNode"]:
        node: Optional["PhyloNode"] = self
        while node is not None:
            yield node
            node = node.parent

    # -- content -----------------------------------------------------------

    def tips(self) -> list["PhyloNode"]:
        if self._tips is None:
            if self.is_tip():
                self._tips = [self]
            else:
                # iterative left-to-right leaf collection (deep-tree safe)
                stack = list(self.children)
                leaves: list["PhyloNode"] = []
                while stack:
                    node = stack.pop()
                    if node.is_tip():
                        leaves.append(node)
                    else:
                        stack.extend(node.children)
                leaves.reverse()
                self._tips = leaves
        return self._tips

    def tip_labels(self) -> list[str]:
        return [t.label for t in self.tips()]

    @property
    def clade(self) -> frozenset[str]:
        """The set of tip labels descending from this node (cached)."""
        if self._clade is None:
            self._clade = frozenset(self.tip_labels())
        return self._clade

    def n_tips(self) -> int:
        return len(self.tips())

    def n_internal(self) -> int:
        return sum(1 for n in self.walk() if not n.is_tip())

    def find_tip(self, name: str) -> Optional["PhyloNode"]:
        for tip in self.tips():
            if tip.label == name:
                return tip
        return None

    def find_internal_by_label(self, label: str) -> Optional["PhyloNode"]:
        for node in self.walk():
            if not node.is_tip() and node.label == label:
                return node
        return None

    def mrca(self, names: Iterable[str]) -> "PhyloNode":
        """Most recent common ancestor of the given tip labels."""
        wanted = set(names)
        missing = wanted - set(self.tip_labels())
        if missing:
            raise KeyError(
                f"taxon not found in tree: {sorted(missing)!r}; "
                f"available tips: {len(self.tip_labels())}"
            )
        best = self
        for node in self.iter_postorder():
            if wanted <= node.clade and node.n_tips() < best.n_tips():
                best = node
        return best

    # -- ids ---------------------------------------------------------------

    def assign_ids(self) -> None:
        """Number tips 1..n then internal nodes n+1.. in pre-order (cladewise).

        The root receives ``n + 1`` and the deepest internal node
        ``n + Nnode`` — the same numbering as ``ape``'s cladewise order and
        the node tables emitted by R3F/MEGA, so ``NodeId``/``Des1``/``Des2``
        are interchangeable across tools.  (Ancestors therefore always have
        *smaller* ids than their descendants; R3F's ``max(Ancestors)`` idiom
        selects the nearest ancestor under this contract.)
        """
        counter = 0
        for node in self.tips():
            counter += 1
            node.node_id = counter
        for node in self.iter_preorder():
            if not node.is_tip():
                counter += 1
                node.node_id = counter

    # -- copies ------------------------------------------------------------

    def copy(self) -> "PhyloNode":
        clone = PhyloNode(
            label=self.label,
            blen=self.blen,
            support=self.support,
            rate=self.rate,
            internal=self.internal,
        )
        clone.node_id = self.node_id
        stack = [(self, clone)]
        while stack:
            src, dst = stack.pop()
            for child in src.children:
                child_clone = PhyloNode(
                    label=child.label,
                    blen=child.blen,
                    support=child.support,
                    rate=child.rate,
                    internal=child.internal,
                )
                child_clone.node_id = child.node_id
                dst.add_child(child_clone, blen=child.blen)
                stack.append((child, child_clone))
        return clone

    def __repr__(self) -> str:  # pragma: no cover - debugging helper
        kind = "tip" if self.is_tip() else "internal"
        return f"PhyloNode({kind}, label={self.label!r}, blen={self.blen})"


def prune_childless_internal_nodes(tree: PhyloNode) -> None:
    """Remove childless nodes that stand for internal (clade) structure.

    After outgroup tips are removed, the internal nodes that anchored them
    have no children left and are pruned bottom-up.  A label does *not*
    spare such a node: the name described the subtree that was just deleted,
    so keeping the node would turn it into a "ghost tip" — ``is_tip()`` is
    ``not children``, so the ingroup would silently gain a taxon.  Genuine
    tips (nodes that were never internal) are untouched.
    """
    changed = True
    while changed:
        changed = False
        for node in list(tree.walk()):
            if node.is_root() or node.children or not node.internal:
                continue
            parent = node.parent
            assert parent is not None
            parent.children.remove(node)
            node.parent = None
            changed = True


def invalidate_caches(root: PhyloNode) -> None:
    """Clear the clade/tip caches of *every* node reachable from *root*.

    Must be called once after any structural mutation that removes or
    re-parents nodes (pruning, contraction, rooting).  Merely clearing the
    mutation point's ancestor chain is not enough because descendant caches
    below other mutated parents would go stale.
    """
    node: Optional[PhyloNode] = root
    while node.parent is not None:
        node = node.parent
    stack = [node]
    while stack:
        n = stack.pop()
        n._clade = None
        n._tips = None
        stack.extend(n.children)


def contract_unary_nodes(
    tree: PhyloNode, discard_root_edge: bool = True
) -> PhyloNode:
    """Contract unary internal nodes, accumulating branch lengths.

    Parameters
    ----------
    discard_root_edge:
        When the root ends up with a single child, replace the root by that
        child and *discard* the connecting edge length.  This mirrors the
        observed behaviour of ``ape::drop.tip`` on outgroup-rooted trees
        (documented in docs/adr/ADR-004-tree-operations.md).

    Returns
    -------
    PhyloNode
        The (possibly new) root of the tree.
    """
    prune_childless_internal_nodes(tree)
    # non-root unary chains
    changed = True
    while changed:
        changed = False
        for node in list(tree.walk()):
            if node.is_root() or len(node.children) != 1:
                continue
            child = node.children[0]
            parent = node.parent
            assert parent is not None
            child.blen = (node.blen or 0.0) + (child.blen or 0.0)
            parent.children.remove(node)
            parent.children.append(child)
            child.parent = parent
            node.parent = None
            node.children = []
            changed = True
    prune_childless_internal_nodes(tree)
    invalidate_caches(tree)
    root = tree
    while discard_root_edge and len(root.children) == 1:
        child = root.children[0]
        child.blen = None
        child.parent = None
        root.children = []
        root = child
        prune_childless_internal_nodes(root)
    invalidate_caches(root)
    return root
