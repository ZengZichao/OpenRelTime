"""Relative lineage rates (:func:`openreltime.rrf_rates` equivalent)."""

from __future__ import annotations

import logging
from typing import Any, Optional

from . import _constants as C
from .rrf import RrfEngine
from .report import RateResult
from .tree import PhyloNode

logger = logging.getLogger("openreltime")


def rrf_rates(
    tree: PhyloNode,
    *,
    mean: str = "geometric",
    rate_ratio_threshold: Optional[float] = C.RATE_RATIO_GUARD_DISABLED,
    **params: Any,
) -> RateResult:
    """Estimate relative lineage rates from a rooted branch-length tree.

    Parameters
    ----------
    tree:
        Rooted binary ingroup tree (as returned by :func:`openreltime.read_tree`).
    mean:
        ``"geometric"`` (default, R3F/msy044 eqs 28-42) or ``"arithmetic"``
        (PNAS 2012 semantics, msy044 eqs 1-27, for cross-validation).
    rate_ratio_threshold:
        Ignored by this function and kept only for signature symmetry with
        :func:`openreltime.rrf_times`.  The R3F reference applies the
        rate-ratio guard in the time pass alone (``rrf_times.R`` lines
        353-373), where it rewrites ``t7.adjust``; ``rrf_rates.R`` has no
        guard and no node times to rewrite, so no threshold can change the
        rates.  Passing a number logs a warning; a warning is also recorded in
        :attr:`RateResult.warnings`.

    Returns
    -------
    RateResult
        Lineage rates per node (tips included; root = 1.0).
    """
    engine = RrfEngine(
        tree, mean=mean, rate_ratio_threshold=rate_ratio_threshold, compute_times=False
    )
    engine.run()

    rates: dict[int, float] = {tree.node_id: 1.0}
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

    warning = ""
    if rate_ratio_threshold is not None:
        warning = (
            "rate_ratio_threshold is ignored by rrf_rates: the R3F rate-ratio "
            "guard (rrf_times.R lines 353-373) only rewrites node times and "
            f"rrf_rates.R applies no guard (threshold {rate_ratio_threshold} "
            "had no effect); use openreltime.rrf_times to apply it"
        )
        logger.warning(
            "rate_ratio_threshold=%s is ignored by rrf_rates (R3F applies the "
            "guard to node times only); the rates are unchanged",
            rate_ratio_threshold,
        )

    merged = {"mean": mean, "rate_ratio_threshold": rate_ratio_threshold, **params}
    return RateResult(
        tree=tree, rates=rates, params=merged, warnings=[warning] if warning else []
    )
