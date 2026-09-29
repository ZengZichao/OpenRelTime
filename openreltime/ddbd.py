"""ddBD: birth-death speciation tree prior from RRF relative times.

Implements the density-based birth-death prior estimator of

    Tao Q, Tamura K, Costa DP, Kumar S.  *A method for assessing
    molecular-clock time calibration.*  ... Bioinformatics 37:i102-i110
    (2021), doi:10.1093/bioinformatics/btab307

following the GPL-3 R3F implementation (``R3F/R/ddbd.R``): the same BD
density, the same initial-value grid (birth/death/sampling), the same SSE and
KL (FNN k-NN, k=5) selection scores and the same L-BFGS-B bounds.  R's
``density()`` kernel estimator and ``FNN::KL.divergence`` are re-implemented
in :func:`_r_density` and :func:`_kl_divergence` so that starting values are
selected identically.
"""

from __future__ import annotations

import logging
import math
from typing import Any, Literal, Optional

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

from . import _constants as C
from .report import DDBDResult
from .tree import PhyloNode
from .times import rrf_times

logger = logging.getLogger("openreltime")

__all__ = ["ddbd"]


def _r_density(x: np.ndarray, n: int = 512, cut: float = 3.0):
    """Gaussian kernel density estimate matching R's ``stats::density``.

    Uses bw.nrd0 (Silverman) and R's default grid ``[min-cut*bw, max+cut*bw]``
    with ``n`` points.
    """
    x = np.asarray(x, dtype=float)
    n_x = x.size
    sd = float(np.std(x, ddof=1))
    q75, q25 = np.quantile(x, [0.75, 0.25])
    iqr = float(q75 - q25)
    lo = min(sd, iqr / 1.349)
    if lo <= 0:
        lo = sd if sd > 0 else (abs(float(x[0])) if x[0] != 0 else 1.0)
    if n_x < 4:
        bw = 0.9 * lo if n_x > 1 else 1.0
    else:
        bw = 0.9 * lo * n_x ** (-0.2)
    grid = np.linspace(float(x.min()) - cut * bw, float(x.max()) + cut * bw, n)
    z = (grid[:, None] - x[None, :]) / bw
    y = norm.pdf(z).sum(axis=1) / (n_x * bw)
    return grid, y


def _kl_divergence(p: np.ndarray, q: np.ndarray, k: int = 5) -> float:
    """1-D k-NN KL divergence estimate matching ``FNN::KL.divergence(p, q, k)``.

    KL*(P||Q) = mean_i log(nu_k(i)/rho_k(i)) + log(m/(n-1))  (d = 1), where
    rho_k is the k-th neighbour distance within *p* and nu_k the k-th
    neighbour distance from p_i to *q*.
    """
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    n, m = p.size, q.size
    k_eff = min(k, m - 1)
    ps = np.sort(p)
    qs = np.sort(q)

    def kth_within(arr: np.ndarray, k: int) -> np.ndarray:
        # k-th nearest neighbour distance among the sample itself (knnDist
        # includes the point's own distance only if duplicated values exist,
        # matching FNN behaviour for tied data)
        out = np.empty(arr.size)
        for i, v in enumerate(arr):
            d = np.abs(arr - v)
            d_sorted = np.sort(d)
            out[i] = d_sorted[min(k, d_sorted.size - 1)]
        return out

    def kth_cross(a: np.ndarray, b: np.ndarray, k: int) -> np.ndarray:
        out = np.empty(a.size)
        for i, v in enumerate(a):
            pos = np.searchsorted(b, v)
            cand: list[float] = []
            for w in b[max(0, pos - k - 1) : pos + k + 1]:
                cand.append(abs(v - w))
            cand.sort()
            out[i] = cand[min(k, len(cand) - 1)]
        return out

    rho = kth_within(ps, k)
    nu = kth_cross(ps, qs, k_eff)
    kl = float(np.mean(np.log(nu / rho)) + math.log(m / (n - 1)))
    return max(kl, 0.0)


def _bd_density(
    t: np.ndarray, birth: float, death: float, rho: float, root_age: float = 1.0
) -> np.ndarray:
    """Birth-death density of relative branching times (ddbd.R lines 414-424)."""
    lam, mu = birth, death
    a = np.exp((mu - lam) * t)
    a1 = math.exp((mu - lam) * root_age)
    prob_t = (rho * (lam - mu)) / (rho * lam + (lam * (1 - rho) - mu) * a)
    prob_t1 = (rho * (lam - mu)) / (rho * lam + (lam * (1 - rho) - mu) * a1)
    vt1 = 1 - (1 / rho) * prob_t1 * a1
    p1 = (1 / rho) * prob_t**2 * a
    return lam * p1 / vt1


def _nll(params: np.ndarray, rel_time: np.ndarray) -> float:
    lam, mu, rho = params
    if lam <= 0 or mu < 0 or rho <= 0 or rho > 1:
        return 1e12
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        gt = _bd_density(rel_time, lam, mu, rho, root_age=1.0)
    if not np.all(np.isfinite(gt)) or np.any(gt <= 0):
        return 1e12
    return float(-np.sum(np.log(gt)))


def ddbd(
    tree: PhyloNode,
    *,
    sampling_frac: Optional[float] = None,
    anchor_node: Optional[int] = None,
    anchor_time: float = 1.0,
    measure: Literal["SSE", "KL"] = "SSE",
    **params: Any,
) -> DDBDResult:
    """Fit a birth-death speciation prior to RRF relative node times.

    Parameters
    ----------
    tree:
        Rooted binary ingroup tree.
    sampling_frac:
        Sampling fraction to *report*.  As in R3F (``ddbd.R`` lines
        491-505) the likelihood optimisation always leaves rho free
        (bounds ``(0, 1)``) and the supplied value only replaces the
        reported sampling fraction; birth and death rates are the optimum of
        the free-rho fit.  Must lie in ``(0, 1]``.
    anchor_node, anchor_time:
        Node id and its absolute time used to rescale relative times.  Keep
        the maximum anchor time <= 10 for numerical stability (R3F guidance).
    measure:
        Initial-value selection score: ``"SSE"`` or ``"KL"``.

    Returns
    -------
    DDBDResult
        Birth rate, death rate and sampling fraction (rates are per the
        anchor time unit).
    """
    if anchor_node is not None and anchor_node <= 0:
        raise ValueError(
            f"anchor node id must be a positive internal-node id, got {anchor_node}"
        )
    times_res = rrf_times(tree, rate_ratio_threshold=C.RATE_RATIO_THRESHOLD)
    internal_times = np.array(
        [times_res.times[n.node_id] for n in tree.walk() if not n.is_tip()]
    )
    rel_time = internal_times / internal_times.max()
    rel_time[rel_time < 0] = 0.0

    if anchor_node is not None:
        ids = [n.node_id for n in tree.walk() if not n.is_tip()]
        if anchor_node not in ids:
            raise KeyError(f"anchor node id {anchor_node} not found; ids: {ids[:5]}...")
        if rel_time[ids.index(anchor_node)] <= 0:
            raise ValueError(
                f"anchor node {anchor_node} has zero relative time and cannot "
                f"be used for rescaling"
            )
        sf = anchor_time / rel_time[ids.index(anchor_node)]
    else:
        sf = 1.0

    if sampling_frac is not None:
        rho0 = float(sampling_frac)
        if not 0.0 < rho0 <= 1.0:
            raise ValueError(
                f"sampling_frac must lie in (0, 1], got {sampling_frac!r}"
            )
    # R3F optimises (lambda, mu, rho) with rho free in both branches
    # (``mle(..., lower=c(0,0,0), upper=c(Inf,Inf,1))``, ddbd.R:469-505); a
    # user-supplied sampling fraction is only substituted when reporting.
    opt_bounds: list[tuple[float, float]] = [(0.0, np.inf), (0.0, np.inf), (0.0, 1.0)]

    # -- initial-value grid (ddbd.R lines 399-403) ---------------------------
    # ``expand.grid(b.rate.try, d.rate.try, s.fr.try)`` varies Var1 (birth)
    # fastest, then death, then sampling; the loop order below reproduces
    # that, so ties in the start-point scores break the same way R3F does.
    grid = [
        (b, d, s)
        for s in C.DDBD_SAMPLING_GRID
        for d in C.DDBD_DEATH_GRID
        for b in C.DDBD_BIRTH_GRID
        if b >= d
    ]
    den_x, den_y = _r_density(rel_time)
    mask = (den_x >= 0) & (den_x <= 1)
    gx, gy = den_x[mask], den_y[mask]

    errs: list[float] = []
    kls: list[float] = []
    for b, d, s in grid:
        bd = _bd_density(gx, b, d, s, root_age=1.0)
        errs.append(float(math.sqrt(np.sum((gy - bd) ** 2))))
        kls.append(_kl_divergence(gy, bd, k=5))

    # stable sort: R3F picks the start with ``match(sorted_score, score)``,
    # i.e. the first grid row among ties
    order = np.argsort(kls if measure == "KL" else errs, kind="stable")

    # -- optimisation: try starts in score order (ddbd.R lines 451-498) ------
    best = None
    for attempt in range(min(C.DDBD_MAX_START_ATTEMPTS, len(order))):
        start = np.array(grid[int(order[attempt])], dtype=float)
        try:
            res = minimize(
                _nll,
                start,
                args=(rel_time,),
                method="L-BFGS-B",
                bounds=opt_bounds,
            )
        except Exception:  # noqa: BLE001 - optimizer guards, mirrors try()
            continue
        # R3F accepts any start that does not raise (``try(...)`` is not a
        # ``try-error``), which includes L-BFGS-B's
        # ABNORMAL_TERMINATION_IN_LNSRCH: the coefficients are usable even
        # though SciPy flags ``success=False``.  Requiring ``res.success`` here
        # would make OpenRelTime give up where R3F returns estimates; only a
        # genuinely non-finite result is rejected.
        if np.all(np.isfinite(res.x)) and np.isfinite(res.fun):
            if not res.success:
                logger.info(
                    "ddBD start %d (%s) terminated abnormally: %s; coefficients "
                    "accepted as R3F does",
                    attempt + 1,
                    tuple(start),
                    res.message,
                )
            best = res.x
            break
    if best is None:
        raise RuntimeError(
            "the best parameter setting cannot be found after trying "
            f"{C.DDBD_MAX_START_ATTEMPTS} starts; check the anchor time unit"
        )

    lam, mu, rho = (float(v) for v in best)
    if sampling_frac is not None:
        b_rate, d_rate, s_frac = lam / sf, mu / sf, rho0
    else:
        b_rate, d_rate, s_frac = lam / sf, mu / sf, rho

    merged = {
        "measure": measure,
        "sampling_frac": sampling_frac,
        "anchor_node": anchor_node,
        # keep the free-rho optimum visible whenever a user-supplied sampling
        # fraction replaced the reported value (R3F semantics)
        "sampling_frac_fitted": rho,
        **params,
    }
    return DDBDResult(
        birth_rate=b_rate,
        death_rate=d_rate,
        sampling_frac=s_frac,
        scale_factor=float(sf),
        anchor_node=anchor_node,
        anchor_time=anchor_time if anchor_node is not None else None,
        params=merged,
    )
