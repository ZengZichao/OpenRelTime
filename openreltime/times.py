"""Relative node times (:func:`openreltime.rrf_times` / ``rrf_rates_times``)."""

from __future__ import annotations

from typing import Any, Optional

from . import _constants as C
from .rrf import RrfEngine
from .report import TimeResult
from .tree import PhyloNode


def rrf_times(
    tree: PhyloNode,
    *,
    mean: str = "geometric",
    rate_ratio_threshold: Optional[float] = C.RATE_RATIO_THRESHOLD,
    normalize: bool = False,
    with_rates: bool = False,
    **params: Any,
) -> TimeResult:
    """Estimate relative node times (and optionally lineage rates).

    Parameters
    ----------
    tree:
        Rooted binary ingroup tree.
    mean:
        ``"geometric"`` (default, R3F-compatible) or ``"arithmetic"``
        (PNAS 2012 semantics; both the local solve and the composite-length
        folding follow the chosen convention).
    rate_ratio_threshold:
        Nodes whose adjusted rate exceeds the threshold (or its reciprocal)
        have their time replaced by the closest clean ancestor's time
        (R3F ``rate.ratio = 20``; the default is
        ``openreltime._constants.RATE_RATIO_THRESHOLD``, the single source of
        registered defaults, cf. docs/parameters.md).  ``None`` disables the
        guard.
    normalize:
        Divide node times by the maximum so the root is 1.0 (MEGA/PNAS 2012
        reporting convention).  Default ``False`` (raw ``t.adj`` like R3F).
    with_rates:
        Also report lineage rates (the ``rrf_rates_times`` mode).

    Returns
    -------
    TimeResult
        Node times (tips = 0) plus rates when ``with_rates`` is set.
    """
    engine = RrfEngine(
        tree, mean=mean, rate_ratio_threshold=rate_ratio_threshold, compute_times=True
    )
    engine.run()

    times: dict[int, float] = {}
    rates: dict[int, float] = {tree.node_id: 1.0}
    for node in tree.walk():
        if node.is_tip():
            times[node.node_id] = 0.0
            continue
        row = engine.row(node)
        times[node.node_id] = float(row.t7a) if row.t7a is not None else 0.0
    for node in tree.walk():
        if node.is_root():
            continue
        parent = node.parent
        assert parent is not None
        prow = engine.row(parent)
        value = prow.r5a if parent.children[0] is node else prow.r6a
        if value is None:
            value = 1.0
        rates[node.node_id] = float(value)

    if normalize and times:
        tmax = max(times.values())
        if tmax > 0:
            times = {nid: t / tmax for nid, t in times.items()}

    merged = {
        "mean": mean,
        "rate_ratio_threshold": rate_ratio_threshold,
        "normalize": normalize,
        "with_rates": with_rates,
        **params,
    }
    return TimeResult(
        tree=tree,
        rates=rates,
        times=times,
        n_rate_guarded=engine.n_exceeded,
        params=merged,
    )


def rrf_rates_times(
    tree: PhyloNode,
    *,
    mean: str = "geometric",
    rate_ratio_threshold: Optional[float] = C.RATE_RATIO_THRESHOLD,
    normalize: bool = False,
    **params: Any,
) -> TimeResult:
    """Estimate relative lineage rates and node times in one pass."""
    return rrf_times(
        tree,
        mean=mean,
        rate_ratio_threshold=rate_ratio_threshold,
        normalize=normalize,
        with_rates=True,
        **params,
    )
