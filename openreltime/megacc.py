"""MEGA-CC (``megacc``) bridge: interoperability layer, optional extra.

This module is *not* on the core algorithmic path.  It provides:

1. a generator for the ``reltimeFromBranchLengths.mao`` analysis-options file
   (INI-style plain text; key semantics per MEGA's ``mparse_mao_file.pas``,
   format verified against the ten real ``.mao`` files shipped with
   RelTime-JA);
2. a runner that invokes ``megacc`` through
   :func:`openreltime._external.run_tool` (validated argument list, no
   shell) and probes its version banner;
3. a parser that reads *megacc's own numbers* — the node ages of its time
   tree, or its tabulated node-time output — back into an R3F-style node
   table whose ``source`` column says where every figure came from.

Provenance contract: :func:`parse_megacc_output` never
runs OpenRelTime's estimator in place of MEGA.  What it cannot read reliably
it refuses, raising :class:`MegaccParseError`, so that an
"OpenRelTime versus MEGA" comparison can never degenerate into a comparison
of OpenRelTime with itself.
"""

from __future__ import annotations

import logging
import re
import tempfile
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from ._external import resolve_executable, run_tool, validate_tokens
from .report import node_table
from .tree import PhyloNode
from .treeio import _extract_newick_from_nexus, parse_newick, read_tree

logger = logging.getLogger("openreltime")

__all__ = ["write_mao", "run_megacc", "parse_megacc_output", "MegaccParseError"]

_MAO_TEMPLATE = """#MAO VERSION=1
[II]ANALYSIS?f=RelTimeFromBranchLengths[II]
[II]ppRelTimeBLens?f=true[II]
"""


def write_mao(path: str | Path, *, overrides: Optional[dict[str, str]] = None) -> Path:
    """Write a minimal ``reltimeFromBranchLengths`` .mao file.

    ``overrides`` adds ``key = value`` lines to the ``[II]`` section after the
    mandatory ``ppRelTimeBLens=true`` line (keys and values are sanitised via
    the shared external-process blacklist).
    """
    lines = _MAO_TEMPLATE.rstrip("\n").split("\n")
    for key, value in (overrides or {}).items():
        validate_tokens(key, value)
        lines.append(f"[II]{key}?f={value}[II]")
    out = Path(path)
    out.write_text("\n".join(lines) + "\n")
    return out


def _safe_path(path: str | Path) -> str:
    text = str(path)
    validate_tokens(text)
    return text


def run_megacc(
    tree: str | Path,
    *,
    megacc: str,
    calibrations: Optional[str | Path] = None,
    outgroup: Optional[list[str]] = None,
    mao: Optional[str | Path] = None,
    workdir: Optional[str | Path] = None,
    **_: Any,
) -> tuple[Path, list[Path]]:
    """Run ``megacc`` RelTime on a branch-length tree.

    Every path is validated and the external process is started by
    :func:`openreltime._external.run_tool` (argument list, no shell).

    Returns
    -------
    (mao_path, output_files)
        The generated ``.mao`` path and the produced output files.
    """
    executable = resolve_executable(megacc)
    work = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="ort_megacc_"))
    work.mkdir(parents=True, exist_ok=True)
    mao_path = Path(mao) if mao else write_mao(work / "reltimeFromBranchLengths.mao")
    tree_token = _safe_path(tree)
    mao_token = _safe_path(mao_path)
    out_token = _safe_path(work / "megacc_out")

    og_pair: list[str] = []
    if outgroup:
        og_file = work / "outgroup.txt"
        og_file.write_text("\n".join(outgroup) + "\n")
        og_pair = ["-g", _safe_path(og_file)]

    calib_pair: list[str] = []
    if calibrations is not None:
        calib_pair = ["-c", _safe_path(calibrations)]

    probe = run_tool([executable, "--version"], check=False)
    banner = (probe.stdout or probe.stderr or b"").decode(errors="replace").strip()
    if banner:
        logger.info("megacc banner: %s", banner.splitlines()[0])

    run_tool(
        [
            executable,
            "-a",
            mao_token,
            "-t",
            tree_token,
            "-o",
            out_token,
            *og_pair,
            *calib_pair,
        ]
    )
    outputs = sorted(p for p in work.iterdir() if p.is_file() and p.suffix != ".mao")
    return mao_path, outputs


#: Value of the ``source`` column added to every parsed table.
MEGACC_SOURCE = "megacc"


#: Header names MEGA/R3F dialects use for the node identifier and its time.
#: Compared after lower-casing and after dropping a parenthesised unit suffix
#: (megacc writes e.g. ``Time (Myr)``).
_ID_COLUMNS = ("nodeid", "node_id", "node", "nodelabel", "node.label")
_TIME_COLUMNS = ("time", "node.time", "age", "divergencetime", "divergence.time", "mean")
#: Optional textual node label column (kept distinct from the numeric id).
_LABEL_COLUMNS = ("nodelabel", "node.label", "label", "treelabel", "nodename")

_UNIT_SUFFIX = re.compile(r"\s*\(.*?\)\s*")


def _normalise_header(column: Any) -> str:
    """``'Time (Myr)'`` -> ``'time'``: lower-case, unit suffix and quotes off."""
    return _UNIT_SUFFIX.sub("", str(column)).strip().strip("\"'").lower()


def _header_map(columns: Any) -> dict[str, str]:
    """Normalised header -> original column name (first occurrence wins)."""
    out: dict[str, str] = {}
    for column in columns:
        out.setdefault(_normalise_header(column), str(column))
    return out


def _numeric_column(raw: pd.DataFrame, names: tuple[str, ...], headers: dict[str, str]):
    """First of *names* that pandas can read as numbers, else ``None``."""
    for name in names:
        column = headers.get(name)
        if column is None:
            continue
        values = pd.to_numeric(raw[column], errors="coerce")
        if int(values.notna().sum()) > 0:
            return column, values
    return None, None


class MegaccParseError(ValueError):
    """megacc left no output whose numbers can be read unambiguously.

    Subclass of :class:`ValueError`.  Raised instead of falling back to
    OpenRelTime's own estimates: the numbers returned by this bridge are
    MEGA's or nothing.
    """


def _tree_from_text(text: str) -> Optional[PhyloNode]:
    """Parse *text* as a Newick/NEXUS tree, or return ``None``.

    A file counts as a time tree only if it actually parses as one (bracket-
    and quote-aware).  Sniffing for ``(``, ``)`` and ``;`` is not enough: a
    megacc log line or a report with parenthesised column headers contains all
    three, so a non-tree would be mistaken for one.
    """
    body = text.strip()
    if not body:
        return None
    if re.search(r"(?im)^\s*TREE\s+[^=]*=", body):
        try:
            body = _extract_newick_from_nexus(body)
        except ValueError:
            return None
    if not body.endswith(";"):
        body += ";"
    try:
        return parse_newick(body)
    except ValueError:  # NewickError is a ValueError subclass
        return None


def _table_delimiter(text: str) -> Optional[str]:
    """Return the delimiter if *text* starts with a node-id + time header."""
    header = next((line for line in text.splitlines() if line.strip()), "")
    if "\t" in header:
        delim = "\t"
    elif "," in header:
        delim = ","
    else:
        return None
    headers = _header_map(header.split(delim))
    if len(headers) < 2:
        return None
    has_id = any(name in _ID_COLUMNS for name in headers)
    has_time = any(name in _TIME_COLUMNS for name in headers)
    return delim if (has_id and has_time) else None


def _node_ages(tree: PhyloNode) -> dict[int, float]:
    """Node ages of a *time tree*, read off megacc's own branch lengths.

    A megacc RelTime tree scales its edges in time, not in substitutions per
    site, so the age of a node is the path length from that node *to the
    present*: the deepest root-to-tip path minus the root-to-node path.  Tips
    therefore come out at 0, matching the ``Time`` column of R3F/MEGA node
    tables.  These are megacc's numbers; nothing is re-estimated here.
    """
    depth: dict[int, float] = {tree.node_id: 0.0}
    stack: list[PhyloNode] = [tree]
    while stack:
        node = stack.pop()
        acc = depth[node.node_id]
        for child in node.children:
            depth[child.node_id] = acc + float(child.blen or 0.0)
            stack.append(child)
    horizon = max(depth.values())
    return {nid: horizon - d for nid, d in depth.items()}


def _frame_from_tree(path: Path, outgroup: Optional[list[str]]) -> pd.DataFrame:
    text = path.read_text()
    fmt = "nexus" if text.lstrip().upper().startswith("#NEXUS") else "newick"
    # read_tree validates (binary, >=3 tips) and returns ape-cladewise ids, so
    # the table is directly comparable with rrf_times(...).to_pandas()
    tree = read_tree(path, fmt=fmt, outgroup=outgroup)
    ages = _node_ages(tree)
    frame = node_table(tree)
    frame.insert(4, "Time", [ages[int(nid)] for nid in frame["NodeId"]])
    frame.insert(5, "source", MEGACC_SOURCE)
    return frame


def _frame_from_table(path: Path, delim: str) -> pd.DataFrame:
    raw = pd.read_csv(path, sep=delim, dtype=str)
    headers = _header_map(raw.columns)
    id_col, node_ids = _numeric_column(raw, _ID_COLUMNS, headers)
    time_col, times = _numeric_column(raw, _TIME_COLUMNS, headers)
    if id_col is None or time_col is None:
        raise MegaccParseError(
            f"{path.name}: cannot identify numeric node-id and time columns "
            f"among {list(raw.columns)!r}; refusing to guess"
        )
    label_col = next(
        (headers[n] for n in _LABEL_COLUMNS if n in headers and headers[n] != id_col),
        None,
    )
    logger.info(
        "megacc table %s read with megacc's own numbering (%s, %s); its node "
        "ids need not match OpenRelTime ids",
        path.name,
        id_col,
        time_col,
    )
    return pd.DataFrame(
        {
            "NodeLabel": raw[label_col].tolist() if label_col else raw[id_col].tolist(),
            "NodeId": node_ids.astype("Int64").tolist(),
            "Time": times.tolist(),
            "source": [MEGACC_SOURCE] * len(raw),
        }
    )


def parse_megacc_output(
    outputs: list[str | Path], outgroup: Optional[list[str]] = None
) -> pd.DataFrame:
    """Parse megacc's own RelTime output into an R3F-style node-time table.

    Parameters
    ----------
    outputs:
        Files produced by :func:`run_megacc` (or any megacc output set).  A
        time tree (Newick/NEXUS, branch lengths = times) is preferred, since
        its node ids are then re-derived with OpenRelTime's ape-cladewise
        numbering and stay comparable with ``rrf_times``; a tabulated node
        table is the fallback and keeps megacc's own numbering (no
        ``Des1``/``Des2``).
    outgroup:
        Tips to remove from the time tree before reporting, so that the
        table covers the ingroup only (as ``rrf_times`` output does).

    Returns
    -------
    pandas.DataFrame
        Always carrying a ``source`` column equal to :data:`MEGACC_SOURCE`:
        every number in it was written by megacc.

    Raises
    ------
    MegaccParseError
        No candidate file exists, more than one candidate is a genuine tree or
        table, or a candidate cannot be turned into a node table.  OpenRelTime's
        estimates are never substituted for megacc's.
    """
    trees: list[Path] = []
    tables: list[tuple[Path, str]] = []
    seen: list[str] = []
    for item in outputs:
        path = Path(item)
        if not path.exists():
            continue
        seen.append(str(path))
        try:
            text = path.read_text(errors="replace")
        except OSError:  # pragma: no cover - unreadable file
            continue
        delim = _table_delimiter(text)
        if delim is not None:
            tables.append((path, delim))
        elif _tree_from_text(text) is not None:
            trees.append(path)
    if not trees and not tables:
        raise MegaccParseError(
            "no usable megacc output found (expected a RelTime time tree in "
            f"Newick/NEXUS form or a tabulated node-time table); inspected: {seen!r}"
        )
    if len(trees) > 1:
        # several genuine trees is an ambiguous input, so stop rather than let
        # the last matching file silently win

        raise MegaccParseError(
            f"several megacc outputs parse as phylogenetic trees: "
            f"{[str(p) for p in trees]!r}; pass exactly one of them"
        )
    if trees:
        try:
            return _frame_from_tree(trees[0], outgroup)
        except MegaccParseError:
            raise
        except Exception as exc:  # noqa: BLE001 - one documented error type
            raise MegaccParseError(
                f"cannot use megacc time tree {str(trees[0])!r}: {exc}"
            ) from exc
    if len(tables) > 1:
        raise MegaccParseError(
            f"several megacc outputs look like node-time tables: "
            f"{[str(p) for p, _ in tables]!r}; pass exactly one of them"
        )
    path, delim = tables[0]
    try:
        return _frame_from_table(path, delim)
    except MegaccParseError:
        raise
    except Exception as exc:  # noqa: BLE001 - turn reader noise into one error
        raise MegaccParseError(
            f"cannot parse megacc table {str(path)!r}: {exc}"
        ) from exc
