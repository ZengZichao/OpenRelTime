"""Tip classification and monophyly assessment.

Provides the taxonomy layer used by outgroup rooting and calibration:

* **User-supplied classification table** (:func:`parse_taxon_table`): a TSV or
  CSV file mapping each tip label to a higher-level group name.
* **Automatic recognition from tip labels**
  (:func:`detect_groups_from_labels`): the convention ``Genus_species`` (also
  ``Genus_species|accession`` etc.) is recognised by splitting each tip label
  at the first separator (``_``, space, ``|`` or ``@``); the prefix becomes
  the group name.
* **Monophyly testing** (:func:`check_monophyly`): a set of tips is
  monophyletic when their MRCA descends to exactly these tips -- the same
  "complete clade" criterion applied to multi-tip outgroups.
* **Token resolution** (:func:`resolve_taxa_tokens`): outgroup and
  calibration tokens may be a tip label *or* a group name; group names are
  expanded to their member tips and warned about when non-monophyletic.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence

import pandas as pd

from .tree import PhyloNode

logger = logging.getLogger("openreltime")

__all__ = [
    "MonophylyResult",
    "check_monophyly",
    "warn_if_not_monophyletic",
    "parse_taxon_table",
    "detect_groups_from_labels",
    "resolve_taxa_tokens",
]

#: Separators recognised by the automatic label-to-group rule.  The genus is
#: the text before the first of these characters (``Genus_species`` is the
#: overwhelming convention in Newick tip labels).
_DEFAULT_SEP = r"[_|@\s]"

# Header aliases accepted by parse_taxon_table (lower-cased).
_TIP_COLUMNS = {"tip", "tips", "label", "tip_label", "name", "species", "taxon_label"}
_GROUP_COLUMNS = {
    "group",
    "groups",
    "clade",
    "taxa",
    "taxon",
    "taxon_group",
    "classification",
    "family",
    "genus",
    "order",
    "class",
}


@dataclass(frozen=True)
class MonophylyResult:
    """Outcome of a monophyly test for one named set of tips.

    ``is_monophyletic`` is true when every requested tip exists in the tree
    *and* the MRCA of those tips descends to no other tip.
    """

    name: str
    tips: tuple[str, ...]
    missing: tuple[str, ...]
    is_monophyletic: bool
    mrca_tips: tuple[str, ...]
    extra_tips: tuple[str, ...]


def check_monophyly(
    tree: PhyloNode, names: Iterable[str], name: str = "group"
) -> MonophylyResult:
    """Test whether *names* form a monophyletic group (complete clade) on *tree*.

    The criterion is exactly the outgroup completeness rule: tracing the
    requested tips back to their most recent common ancestor, *all* tips
    descending from that ancestor must belong to the requested set.
    """
    labels = set(tree.tip_labels())
    requested = list(dict.fromkeys(str(n) for n in names))
    found = [n for n in requested if n in labels]
    missing = [n for n in requested if n not in labels]
    if not found:
        return MonophylyResult(
            name=name,
            tips=(),
            missing=tuple(missing),
            is_monophyletic=False,
            mrca_tips=(),
            extra_tips=(),
        )
    mrca = tree.mrca(found)
    mrca_tips = tuple(sorted(mrca.tip_labels()))
    extra = sorted(set(mrca_tips) - set(found))
    return MonophylyResult(
        name=name,
        tips=tuple(found),
        missing=tuple(missing),
        is_monophyletic=not missing and not extra,
        mrca_tips=mrca_tips,
        extra_tips=tuple(extra),
    )


def warn_if_not_monophyletic(
    tree: PhyloNode, names: Iterable[str], name: str = "group"
) -> MonophylyResult:
    """Run :func:`check_monophyly` and log a warning when it fails."""
    result = check_monophyly(tree, names, name=name)
    if not result.is_monophyletic:
        if result.missing:
            logger.warning(
                "taxon group %r contains names absent from the tree: %s",
                name,
                ", ".join(result.missing),
            )
        if result.extra_tips:
            shown = ", ".join(result.extra_tips[:8])
            if len(result.extra_tips) > 8:
                shown += f", ... (+{len(result.extra_tips) - 8})"
            logger.warning(
                "taxon group %r is NOT monophyletic: the MRCA of its %d tips "
                "also descends to %d additional tip(s) (%s), so the group is "
                "not a complete clade",
                name,
                len(result.tips),
                len(result.extra_tips),
                shown,
            )
    return result


def parse_taxon_table(path: str | Path) -> dict[str, list[str]]:
    """Read a tip-classification table mapping group names to tip labels.

    The file is a TSV (or CSV) with a header row containing one tip column
    (``tip``/``label``/``name``/... ) and one group column
    (``group``/``clade``/``taxa``/``classification``/...).  With exactly two
    columns and unrecognised headers, the first column is taken as the tip
    label and the second as the group name.  A tip may appear once; a group
    collects all of its tips.  Returns ``{group_name: [tip, ...]}``.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"taxon table not found: {path}")
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else None
    df = pd.read_csv(
        path, sep=sep, dtype=str, keep_default_na=False, engine="python"
    )
    df = df.rename(columns={c: str(c).strip().lower() for c in df.columns})
    tip_col = next((c for c in df.columns if c in _TIP_COLUMNS), None)
    group_col = next((c for c in df.columns if c in _GROUP_COLUMNS), None)
    if tip_col is None or group_col is None:
        if len(df.columns) == 2:
            tip_col, group_col = df.columns[0], df.columns[1]
        else:
            raise ValueError(
                f"taxon table {path} must have a tip column "
                f"(one of {sorted(_TIP_COLUMNS)}) and a group column "
                f"(one of {sorted(_GROUP_COLUMNS)}); found {list(df.columns)}"
            )
    table: dict[str, list[str]] = {}
    seen: dict[str, str] = {}
    for _, row in df.iterrows():
        tip = str(row[tip_col]).strip()
        group = str(row[group_col]).strip()
        if not tip or not group or group in (".", "NA", "-"):
            continue
        if tip in seen:
            raise ValueError(
                f"tip {tip!r} appears twice in taxon table {path} "
                f"(groups {seen[tip]!r} and {group!r})"
            )
        seen[tip] = group
        table.setdefault(group, []).append(tip)
    if not table:
        raise ValueError(f"taxon table {path} contains no usable rows")
    return table


def detect_groups_from_labels(
    labels: Iterable[str],
    *,
    sep: str = _DEFAULT_SEP,
    field: int = 0,
) -> dict[str, list[str]]:
    """Infer groups from tip labels by splitting at the first separator.

    ``Homo_sapiens`` -> group ``Homo``; ``Pan_paniscus|ABC1`` -> ``Pan``.
    A label without any separator forms a singleton group of itself.  Returns
    ``{group_name: [tip, ...]}`` with tips in input order.
    """
    pattern = re.compile(sep)
    table: dict[str, list[str]] = {}
    for label in labels:
        label = str(label)
        parts = pattern.split(label.strip(), maxsplit=1)
        group = parts[field] if len(parts) > field else label
        group = group.strip()
        if not group:
            group = label
        table.setdefault(group, []).append(label)
    return table


def resolve_taxa_tokens(
    tree: PhyloNode,
    tokens: Sequence[str],
    *,
    taxon_table: Optional[Mapping[str, Sequence[str]]] = None,
    context: str = "taxon set",
) -> list[str]:
    """Resolve a list of name tokens to tip labels, expanding group names.

    Each token is resolved as a tip label first; when it is not a tip of the
    tree it is looked up as a *group name* in ``taxon_table`` when given, or
    otherwise against groups auto-detected from the tip labels.  Expanded
    groups are checked for monophyly and a warning is logged when the check
    fails (the caller decides whether that is fatal, e.g. the outgroup
    completeness check in :func:`openreltime.treeio.root_outgroup`).
    """
    labels = set(tree.tip_labels())
    resolved: list[str] = []
    for token in tokens:
        token = str(token).strip()
        if not token:
            continue
        if token in labels:
            resolved.append(token)
            continue
        table = taxon_table
        source = "taxon table"
        if table is None:
            table = detect_groups_from_labels(tree.tip_labels())
            source = "tip-label prefixes (auto-detected)"
        if token in table:
            members = [m for m in table[token] if m in labels]
            absent = sorted(set(table[token]) - set(members))
            if absent:
                logger.warning(
                    "%s %r: %d of %d member tips of group %r are absent from "
                    "the tree: %s",
                    context,
                    token,
                    len(absent),
                    len(table[token]),
                    token,
                    ", ".join(absent),
                )
            if not members:
                raise KeyError(
                    f"{context} token {token!r}: none of its group members "
                    f"were found in the tree"
                )
            logger.info(
                "%s token %r expanded to group with %d tips via %s",
                context,
                token,
                len(members),
                source,
            )
            warn_if_not_monophyletic(tree, members, name=token)
            resolved.extend(m for m in members if m not in resolved)
        else:
            raise KeyError(
                f"{context} taxon {token!r} not found in tree (as tip or "
                f"group); check the spelling or provide a taxon table"
            )
    if not resolved:
        raise KeyError(f"{context} is empty after resolution")
    return resolved
