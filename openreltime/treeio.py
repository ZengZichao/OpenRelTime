"""Tree input/output: Newick & NEXUS parsing, outgroup rooting, validation."""

from __future__ import annotations

import logging
import math
import re
from pathlib import Path
from typing import Optional, Sequence

from .taxonomy import resolve_taxa_tokens
from .tree import PhyloNode, contract_unary_nodes, invalidate_caches

logger = logging.getLogger("openreltime")

# ---------------------------------------------------------------------------
# Newick parsing
# ---------------------------------------------------------------------------


class NewickError(ValueError):
    """Raised when a Newick string cannot be parsed."""


_NEWICK_TOKEN = re.compile(
    r"""
    \s*(?:
        (?P<open>\()
      | (?P<close>\))
      | (?P<comma>,)
      | (?P<semi>;)
      | (?P<label>[^,:;()\[\]]+)
    )
    """,
    re.VERBOSE,
)

_NUMBER = re.compile(r"\s*([-+]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][-+]?[0-9]+)?)")

_QUOTES = "'\""

#: ``[&rate=0.123]``-style values inside a Newick/NEXUS comment.  The key is
#: matched case-insensitively and may sit anywhere in a comma-separated
#: attribute list (FigTree writes ``[&rate=…,height_95%_max=…]``).
_RATE_ATTR = re.compile(r"(?i)\brate\s*=\s*([-+]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][-+]?[0-9]+)?)")


def _read_comment(text: str, pos: int) -> tuple[str, int]:
    """Read a ``[...]`` Newick comment starting at *pos*.

    Bracket-aware (NHOL allows nested comments, e.g. ``[a [b] c]``) and
    quote-aware (``[it's a comment]``), so that unlike a ``\\[[^\\]]*\\]``
    regex it can never cut through a label.  Returns the comment body
    (without the delimiters) and the position after the closing bracket.
    """
    assert text[pos] == "["
    depth = 0
    start = pos + 1
    while pos < len(text):
        ch = text[pos]
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
            if depth == 0:
                return text[start : pos], pos + 1
        elif ch in _QUOTES:
            _, pos = _read_quoted(text, pos)
            continue
        pos += 1
    raise NewickError("unterminated [...] comment in Newick string")


def _comment_rate(body: str) -> Optional[float]:
    """Return the ``rate`` value of an NH annotation comment, if any."""
    match = _RATE_ATTR.search(body)
    return float(match.group(1)) if match else None


def _read_quoted(text: str, pos: int) -> tuple[str, int]:
    """Read a quoted label starting at *pos*; doubled quotes are escapes."""
    quote = text[pos]
    pos += 1
    out: list[str] = []
    while pos < len(text):
        ch = text[pos]
        if ch == quote:
            if pos + 1 < len(text) and text[pos + 1] == quote:
                out.append(quote)
                pos += 2
                continue
            return "".join(out), pos + 1
        out.append(ch)
        pos += 1
    raise NewickError("unterminated quoted label in Newick string")


def parse_newick(text: str) -> PhyloNode:
    """Parse a Newick string into a :class:`~openreltime.tree.PhyloNode`.

    The parser is an explicit-stack state machine (no recursion), so trees
    far deeper than Python's recursion limit parse fine, and quoted labels
    (``'A, B'``; a literal quote is written doubled) may contain Newick
    punctuation, square brackets included.

    ``[...]`` comments (NHOL / FigTree / R3F dialect) are consumed by a
    bracket-aware scanner rather than a regex, so they neither corrupt
    labels nor hide from the parser.  A ``[&rate=…]`` annotation is kept on
    the node as :attr:`~openreltime.tree.PhyloNode.rate`; every other
    attribute (``[&R]``, colours, support comments) is discarded.
    """
    text = text.strip()
    if not text.endswith(";"):
        raise NewickError("Newick string must end with ';'")

    class _Cursor:
        __slots__ = ("pos",)

        def __init__(self) -> None:
            self.pos = 0

        def error(self, message: str) -> NewickError:
            snippet = text[max(0, self.pos - 25) : self.pos + 25]
            return NewickError(f"{message} near position {self.pos}: {snippet!r}")

    cur = _Cursor()

    def skip_ws() -> None:
        while cur.pos < len(text) and text[cur.pos] in " \t\r\n":
            cur.pos += 1

    def peek() -> Optional[str]:
        skip_ws()
        return text[cur.pos] if cur.pos < len(text) else None

    def skip_comments(node: Optional[PhyloNode] = None) -> None:
        """Consume every ``[...]`` comment at the cursor, keeping rates."""
        while True:
            skip_ws()
            if cur.pos >= len(text) or text[cur.pos] != "[":
                return
            body, cur.pos = _read_comment(text, cur.pos)
            if node is None:
                continue
            rate = _comment_rate(body)
            if rate is not None:
                node.rate = rate

    def parse_label_name(fragment: str) -> Optional[str]:
        clean = fragment.strip().strip("'\"")
        return clean or None

    def read_bare_label() -> str:
        m = _NEWICK_TOKEN.match(text, cur.pos)
        if not m or m.group("label") is None:
            raise cur.error("expected a taxon label")
        cur.pos = m.end()
        name = parse_label_name(m.group("label"))
        if name is None:
            raise cur.error("tip without a name")
        return name

    def read_label() -> str:
        if text[cur.pos] in _QUOTES:
            name, cur.pos = _read_quoted(text, cur.pos)
            return name
        return read_bare_label()

    def read_branch_length(node: PhyloNode) -> None:
        skip_ws()
        if cur.pos < len(text) and text[cur.pos] == ":":
            cur.pos += 1
            m2 = _NUMBER.match(text, cur.pos)
            if not m2:
                raise cur.error("invalid branch length")
            cur.pos = m2.end()
            node.blen = float(m2.group(1))

    def read_node_decoration(node: PhyloNode) -> None:
        """Read ``[comment] label :blen [comment] label`` after a ``)``."""
        node.internal = True
        skip_comments(node)
        if cur.pos < len(text) and text[cur.pos] in _QUOTES:
            label, cur.pos = _read_quoted(text, cur.pos)
            if label:
                if _looks_like_number(label):
                    node.support = float(label)
                else:
                    node.label = label
        else:
            m = _NEWICK_TOKEN.match(text, cur.pos)
            if m and m.group("label") is not None:
                cur.pos = m.end()
                label = parse_label_name(m.group("label"))
                if label is not None:
                    if _looks_like_number(label):
                        node.support = float(label)
                    else:
                        node.label = label
        read_branch_length(node)
        skip_comments(node)

    root: Optional[PhyloNode] = None
    stack: list[PhyloNode] = []
    skip_comments()  # leading NHOL comments, e.g. the NEXUS "[&R]" marker
    while True:
        skip_ws()
        if cur.pos >= len(text):
            raise cur.error("unexpected end of Newick string")
        ch = text[cur.pos]

        if ch == "[":
            skip_comments()
            continue

        if ch == "(":
            node = PhyloNode(internal=True)
            if stack:
                stack[-1].add_child(node)
            stack.append(node)
            cur.pos += 1
            continue

        if ch == ")":
            if not stack:
                raise cur.error("unexpected ')'")
            node = stack.pop()
            if not node.children:
                raise cur.error("expected a taxon label")
            cur.pos += 1
            read_node_decoration(node)
            if stack:
                nxt = peek()
                if nxt == ",":
                    cur.pos += 1
                    continue
                if nxt in (")", "["):
                    continue
                raise cur.error("expected ',' or ')' inside a clade")
            root = node
            break

        if ch in ",:":
            raise cur.error("expected a taxon label")
        tip = PhyloNode(label=read_label())
        skip_comments(tip)
        read_branch_length(tip)
        skip_comments(tip)
        if stack:
            stack[-1].add_child(tip)
            nxt = peek()
            if nxt == ",":
                cur.pos += 1
                continue
            if nxt in (")", "["):
                continue
            raise cur.error("expected ',' or ')' inside a clade")
        root = tip
        break

    skip_ws()
    if peek() != ";":
        raise cur.error("expected ';' terminator")
    assert root is not None
    return root


def _looks_like_number(txt: str) -> bool:
    try:
        float(txt)
        return True
    except ValueError:
        return False


# ---------------------------------------------------------------------------
# NEXUS support
# ---------------------------------------------------------------------------

def _extract_newick_from_nexus(text: str) -> str:
    """Pull the last TREE statement out of a NEXUS document.

    Handles ``BEGIN TREES;`` blocks with ``[&R]`` root markers and
    ``[&rate=...]`` annotations (R3F / FigTree dialect).  Annotations are
    *kept*: :func:`parse_newick` consumes the comments itself and restores
    the per-node rates, so a rate written by :func:`write_nexus` survives a
    read-back.
    """
    trees = re.findall(r"(?im)^\s*TREE\s+[^=]*=\s*(.+?)\s*$", text)
    if not trees:
        raise ValueError("no TREE statement found in NEXUS content")
    return trees[-1].strip().rstrip(";")


def node_rates(tree: PhyloNode) -> dict[int, float]:
    """Collect the ``[&rate=...]`` annotations carried by *tree*'s nodes.

    Keys are ``node_id`` (see :meth:`~openreltime.tree.PhyloNode.assign_ids`),
    so the result is directly usable as the ``rate_map`` of
    :func:`write_nexus` / :func:`_emit_newick`.  Nodes without an annotation
    are omitted.
    """
    return {n.node_id: float(n.rate) for n in tree.walk() if n.rate is not None}


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

def _format_number(value: float) -> str:
    if value == 0:
        return "0.0"
    return f"{value:.10g}"


def _escape_label(label: str) -> str:
    if re.search(r"[\s(),:;\[\]']", label):
        return "'" + label.replace("'", "''") + "'"
    return label


def _emit_newick(
    node: PhyloNode,
    *,
    rate_map: Optional[dict[int, float]] = None,
    time_map: Optional[dict[int, float]] = None,
    include_internal_labels: bool = False,
) -> str:
    """Serialise a subtree to Newick text (with trailing semicolon).

    Explicit-stack emitter: safe for trees deeper than the recursion limit.
    ``time_map`` maps ``node_id`` to node *time*; when given, branch lengths
    are emitted as ``time(parent) - time(child)`` (the timetree convention of
    R3F), ignoring ``blen``.  ``rate_map`` adds ``[&rate=...]`` annotations.
    """
    rates = rate_map or {}

    def blen_text(n: PhyloNode) -> str:
        if n.is_root():
            return ""
        assert n.parent is not None
        if time_map is not None:
            length = time_map[n.parent.node_id] - time_map[n.node_id]
        else:
            length = n.blen or 0.0
        return f":{_format_number(length)}"

    def label_text(n: PhyloNode) -> str:
        if n.is_tip():
            return _escape_label(n.label or "")
        if n.label and include_internal_labels:
            # write_nexus passes include_internal_labels=True (it always kept
            # internal labels); plain to_newick only on request
            return _escape_label(n.label)
        return ""

    def annotate(n: PhyloNode, label: str) -> str:
        if n.node_id in rates:
            return f"{label}[&rate={_format_number(rates[n.node_id])}]"
        return label

    parts: list[str] = []
    stack: list[tuple[Optional[PhyloNode], int]] = [(node, 0)]
    while stack:
        n, state = stack.pop()
        if state == 2:  # comma between siblings
            parts.append(",")
            continue
        if n.is_tip():
            parts.append(annotate(n, label_text(n)) + blen_text(n))
            continue
        if state == 0:
            parts.append("(")
            stack.append((n, 1))
            for i, child in enumerate(reversed(n.children)):
                if i > 0:
                    stack.append((None, 2))
                stack.append((child, 0))
        else:  # state == 1: closing an internal node
            parts.append(")")
            parts.append(annotate(n, label_text(n)))
            parts.append(blen_text(n))
    return "".join(parts) + ";"


def to_newick(
    node: PhyloNode,
    *,
    include_internal_labels: bool = False,
    time_map: Optional[dict[int, float]] = None,
) -> str:
    """Serialise a subtree to Newick (with trailing semicolon).

    ``time_map`` maps ``node_id`` to node *time*; when given, branch lengths
    are emitted as ``time(parent) - time(child)`` (the timetree convention of
    R3F), ignoring ``blen``.
    """
    return _emit_newick(node, time_map=time_map, include_internal_labels=include_internal_labels)


def write_nexus(
    tree: PhyloNode,
    *,
    rate_map: Optional[dict[int, float]] = None,
    time_map: Optional[dict[int, float]] = None,
    title: str = "ORT_tree",
) -> str:
    """Serialise a tree to NEXUS with ``[&rate=...]`` annotations.

    ``rate_map`` maps ``node_id`` to a relative rate; tips and internal nodes
    receive ``[&rate=...]`` comments like R3F's ``rrf_rates_times`` output.
    The TAXA block stays annotation-free (R3F strips it there as well).
    """
    newick_body = _emit_newick(
        tree, rate_map=rate_map, time_map=time_map, include_internal_labels=True
    )
    return "\n".join(
        [
            "#NEXUS",
            f"[OpenRelTime {title}]",
            "BEGIN TAXA;",
            f"    DIMENSIONS NTAX = {tree.n_tips()};",
            "    TAXLABELS",
            *[f"        {_escape_label(lbl)}" for lbl in tree.tip_labels()],
            "        ;",
            "END;",
            "BEGIN TREES;",
            f"    TREE {title} = [&R] {newick_body};",
            "END;",
            "",
        ]
    )


# ---------------------------------------------------------------------------
# Reading with validation
# ---------------------------------------------------------------------------

def read_tree(
    path: str | Path,
    fmt: str = "newick",
    outgroup: Optional[Sequence[str]] = None,
    *,
    resolve_polytomy: str = "error",
    seed: Optional[int] = None,
    taxon_table: Optional[dict[str, list[str]]] = None,
    outgroup_check: str = "error",
) -> PhyloNode:
    """Read, validate, optionally root on an outgroup and prune it.

    Parameters
    ----------
    path:
        Tree file.
    fmt:
        ``"newick"`` (default) or ``"nexus"``.
    outgroup:
        Tip names defining the rooting outgroup; the tree is rooted on their
        MRCA and the outgroup tips are removed (output contains ingroup only).
        Each entry may also be a *group name* (resolved through
        ``taxon_table`` or against prefixes auto-detected from the tip
        labels, e.g. ``"Mus"`` for every ``Mus_...`` tip).
    resolve_polytomy:
        ``"error"`` (default) refuses non-binary trees; ``"random"`` resolves
        polytomies with zero-length branches (a warning is logged).
    seed:
        Seed for ``resolve_polytomy="random"``.
    taxon_table:
        Optional ``{group_name: [tip, ...]}`` mapping (see
        :func:`openreltime.taxonomy.parse_taxon_table`) used to expand
        group-name outgroup tokens.
    outgroup_check:
        What to do when a multi-tip outgroup does **not** form a complete
        clade (its MRCA descends to tips outside the outgroup): ``"error"``
        (default) raises :class:`ValueError`; ``"warn"`` logs a warning and
        roots on the MRCA anyway.

    Returns
    -------
    PhyloNode
        Rooted, binary, outgroup-free tree with ids assigned.  Nodes that
        carried a ``[&rate=...]`` annotation keep it as ``node.rate``
        (see :func:`node_rates`).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"tree file not found: {path}")
    text = path.read_text()
    if fmt.lower() == "nexus":
        newick = _extract_newick_from_nexus(text)
    elif fmt.lower() in ("newick",):
        newick = text.strip()
    else:
        raise ValueError(f"unsupported format: {fmt!r} (use 'newick' or 'nexus')")

    # NB: no comment stripping here -- parse_newick consumes ``[...]`` comments
    # itself and keeps ``[&rate=...]`` on the node, which a blanket regex
    # would silently destroy (and which corrupts quoted labels containing ``[``).
    newick = newick.strip()
    if not newick.endswith(";"):
        newick += ";"
    tree = parse_newick(newick)
    _validate(tree, resolve_polytomy, seed)

    if outgroup:
        if taxon_table is not None or any(
            tree.find_tip(str(tok)) is None for tok in outgroup
        ):
            outgroup = resolve_taxa_tokens(
                tree, outgroup, taxon_table=taxon_table, context="outgroup"
            )
        tree = root_outgroup(tree, outgroup, on_incomplete=outgroup_check)
        _validate(tree, resolve_polytomy, seed)
    tree.assign_ids()
    return tree


def _validate(tree: PhyloNode, resolve_polytomy: str, seed: Optional[int]) -> None:
    labels = [t.label for t in tree.tips()]
    duplicates = {lbl for lbl in labels if labels.count(lbl) > 1}
    if duplicates:
        raise ValueError(f"duplicate tip names: {sorted(duplicates)!r}")
    for node in tree.walk():
        if node.blen is None and not node.is_root():
            # missing lengths are treated as zero by convention; make it loud
            logger.warning(
                "branch length missing at node %r; assuming 0.0",
                node.label or node.node_id,
            )
            node.blen = 0.0
        if node.blen is not None and not math.isfinite(node.blen):
            raise ValueError(f"branch length at node {node.label!r} is NaN/Inf")
        if node.blen is not None and node.blen < 0:
            raise ValueError(
                f"negative branch length at node {node.label!r}; check the input tree"
            )
    tips = tree.tips()
    if len(tips) < 3:
        raise ValueError(
            f"tree has {len(tips)} tips; at least 3 are required "
            "(a 2-taxon ingroup has no resolvable rate solution)"
        )
    # A node that was read or built as an internal node must not survive as a
    # leaf: after outgroup removal such a leftover ("ghost tip") would be
    # counted as ingroup membership and named like a taxon.
    ghosts = sorted(
        str(n.label) if n.label else f"node{n.node_id}"
        for n in tips
        if n.internal
    )
    if ghosts:
        raise ValueError(
            f"{len(ghosts)} internal node(s) degenerated into tips "
            f"({', '.join(ghosts[:5])}); an internal node whose descendants "
            "were removed must be pruned, not kept as a taxon "
            "(see openreltime.tree.prune_childless_internal_nodes)"
        )
    if not tree.is_binary():
        if resolve_polytomy == "random":
            _resolve_polytomies_randomly(tree, seed)
        else:
            spots = [
                (n.label or f"node{n.node_id}")
                for n in tree.walk()
                if len(n.children) > 2
            ]
            raise ValueError(
                f"tree contains polytomies at {spots!r}; only binary trees are "
                "allowed -- use resolve_polytomy='random' (--resolve random) "
                "to binarise with zero-length branches"
            )
    n_nodes = sum(1 for _ in tree.walk())
    n_edges = n_nodes - 1
    zero = sum(1 for n in tree.walk() if not n.is_root() and n.blen == 0.0)
    if n_edges and zero / n_edges >= 0.10:
        logger.warning(
            "there are too many zero-length branches (>= %.0f%%) in the tree; "
            "rates and node times on those branches may not be reliable",
            100.0 * zero / n_edges,
        )


def _resolve_polytomies_randomly(tree: PhyloNode, seed: Optional[int]) -> None:
    import random

    # Polysemy must stay reproducible across runs, so a seeded generator is
    # required here; this is not a security context.
    rng = random.Random(seed)
    resolved: list[str] = []
    for node in list(tree.walk()):
        while len(node.children) > 2:
            picked = rng.sample(node.children, 2)
            new_node = PhyloNode(blen=0.0)
            for child in picked:
                node.children.remove(child)
                new_node.add_child(child, blen=child.blen)
            node.add_child(new_node, blen=0.0)
            resolved.append(node.label or f"node{node.node_id}")
    if resolved:
        logger.warning(
            "polytomies resolved randomly with zero-length branches at: %s",
            ", ".join(sorted(resolved)),
        )


# ---------------------------------------------------------------------------
# Outgroup rooting (ape::root(resolve.root=TRUE) + drop.tip equivalent)
# ---------------------------------------------------------------------------

def check_outgroup_clade(tree: PhyloNode, outgroup: Sequence[str]) -> list[str]:
    """Return tips outside *outgroup* that descend from the outgroup MRCA.

    An empty list means the outgroup forms a complete clade: tracing the
    specified branches back to their most recent common ancestor, every tip
    under that ancestor belongs to the outgroup.  A single-tip outgroup is
    always complete.
    """
    names = list(dict.fromkeys(str(n) for n in outgroup))
    if len(names) <= 1:
        return []
    mrca = tree.mrca(names)
    return sorted(set(mrca.tip_labels()) - set(names))


def root_outgroup(
    tree: PhyloNode, outgroup: Sequence[str], *, on_incomplete: str = "error"
) -> PhyloNode:
    """Root *tree* on the MRCA of *outgroup* and remove the outgroup tips.

    Semantics aligned with ``ape::root(..., resolve.root = TRUE)`` followed by
    ``ape::drop.tip`` (verified against ape 5.8; see docs/adr/ADR-004):

    * a new root is inserted with the outgroup MRCA as one child.  The
      outgroup side carries a zero-length stem, except when the MRCA is a
      single tip, which keeps its own terminal branch;
    * edges along the reversed remaining path keep their lengths and unary
      nodes on that path are contracted (accumulating lengths);
    * after outgroup removal, unary roots left behind are contracted and the
      final edge above a single-child root is discarded.

    A multi-tip outgroup must form a complete clade: all tips under its MRCA
    must belong to the outgroup.  ``on_incomplete="error"`` (default) raises
    :class:`ValueError` otherwise; ``"warn"`` logs a warning and roots on the
    MRCA anyway.  The ingroup must retain at least three tips.
    """
    missing = [name for name in outgroup if tree.find_tip(name) is None]
    if missing:
        raise KeyError(
            f"outgroup taxon {missing!r} not found in tree; check the spelling"
        )
    mrca = tree.mrca(outgroup)
    og_set = set(outgroup)

    extra = check_outgroup_clade(tree, outgroup)
    if extra:
        shown = ", ".join(extra[:8]) + (f", ... (+{len(extra) - 8})" if len(extra) > 8 else "")
        message = (
            f"outgroup {sorted(og_set)!r} is not a complete clade: the MRCA of "
            f"its {len(og_set)} tips also descends to {len(extra)} additional "
            f"tip(s) ({shown}); a multi-tip outgroup must contain every tip "
            f"below its MRCA (monophyly/complete-clade requirement)"
        )
        if on_incomplete == "error":
            raise ValueError(message)
        logger.warning("%s; rooting on the MRCA anyway", message)

    if not mrca.is_root():
        new_root = PhyloNode()
        orig_parent = mrca.parent
        mrca_stem = mrca.blen or 0.0  # edge above the outgroup MRCA
        if orig_parent is not None:
            orig_parent.children.remove(mrca)
        # outgroup side: tip keeps its terminal branch; an internal MRCA gets
        # a zero-length stem (ape resolve.root behaviour)
        if mrca.is_tip():
            new_root.add_child(mrca, blen=mrca_stem)
        else:
            new_root.add_child(mrca, blen=0.0)
        # reversed path from the old parent toward the old root; each chain
        # node carries the length of the edge above the previous path node
        prev_path_node: Optional[PhyloNode] = mrca
        prev_chain: Optional[PhyloNode] = None
        node: Optional[PhyloNode] = orig_parent
        while node is not None:
            if prev_path_node is mrca:
                edge = 0.0 if mrca.is_tip() else mrca_stem
            else:
                edge = prev_path_node.blen or 0.0
            chain = PhyloNode(label=node.label, blen=edge)
            for child in list(node.children):
                if child is prev_path_node:
                    continue
                chain.add_child(child, blen=child.blen)
            if prev_chain is None:
                new_root.add_child(chain, blen=edge)
            else:
                prev_chain.add_child(chain, blen=edge)
            prev_chain = chain
            prev_path_node = node
            node = node.parent
        tree = contract_unary_nodes(new_root, discard_root_edge=False)

    # -- drop outgroup tips --------------------------------------------------
    for name in sorted(og_set):
        tip = tree.find_tip(name)
        if tip is None or tip.is_root():
            continue
        parent = tip.parent
        assert parent is not None
        parent.children.remove(tip)
        tip.parent = None
    invalidate_caches(tree)
    tree = contract_unary_nodes(tree, discard_root_edge=True)

    if tree.n_tips() < 3:
        raise ValueError(
            f"ingroup has {tree.n_tips()} tips after outgroup removal; "
            "at least 3 are required (add taxa or adjust the outgroup)"
        )
    return tree
