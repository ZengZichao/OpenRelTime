"""Analytical confidence intervals for calibrated divergence times.

Implements the delta method of

    Tao Q, Tamura K, Mello B, Kumar S.  *Reliable confidence intervals for
    RelTime estimates of evolutionary divergence times.*  Mol Biol Evol
    37:280-290 (2020), doi:10.1093/molbev/msz236

* branch variance decomposition ``v(b) = vS(b) + vR(b)`` (eq 8);
* ``RV(R)`` recovered from the observed rate variance minus the sampling
  contribution (eqs 9-13), with derivatives of the RRF rate function (eqs 1-4)
  taken by central differences (step ``h = max(1e-8, 1e-5*|b|)``);
* rate-heterogeneity variance per branch ``vR(b_j) = b_j^2 / sum_j(b_j^2) *
  RV(R)`` (eq 14), so that ``sum_j vR(b_j) = RV(R)``;
* node-time variance ``v(t_i) = sum_j (dt_i/db_j)^2 v(b_j)`` (eq 7).  The
  paper evaluates this with closed forms for shallow nodes and a recursion
  for deeper ones (eqs 15-17); here the derivatives of the *whole-pipeline*
  time function are taken numerically, which is algebraically identical to
  the recursion under the no-covariance assumption and stays valid for
  arbitrary calibrated paths (see docs/methods.md and ADR-005, indexed in docs/adr/README.md).
* CIs use the final (calibrated) rates and are truncated at the imposed
  hard calibration bounds.

Two scale/consistency rules are enforced because the printed equations mix
them:

1. *Per-lineage scale.*  eq 9 averages over lineages while eqs 11-13 sum
   ``sv(r_j)`` over them, so subtracting the raw sum from ``Vobs(R)`` clamps
   ``RV(R)`` to 0 on every non-trivial tree and the whole rate-heterogeneity
   component disappears.  Both terms are therefore kept on the per-lineage
   scale of eq 9: ``RV(R) = Vobs(R) - SV(R)/N``.
2. *One estimator throughout.*  The differentiation pipeline re-runs the
   averaging convention (``mean``), the rate-ratio guard and the calibration
   bounds *actually used* for the point estimate, so a ``mean="arithmetic"`` or
   ``method="effective"`` result is qualified by derivatives of that same
   estimator.

A calibrated node's age is a fixed point of its own calibration - either the
hard bound itself (an attractor of the violation loop, so t(b+h) = t(b-h) =
bound) or the midpoint of the feasible interval (t = (min+max)/2 whatever the
branch lengths are).  Its derivative through eq (7) is therefore exactly zero,
which would report "zero uncertainty in age".  Such a node is instead given the
uncertainty its calibration asserts (effective bounds, density quantiles, or
the user-supplied min-max pair), and every row whose interval does not come
from eq (7) is flagged in the ``se_reliable`` column with an explanatory
``notes`` entry.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Callable, Optional

import numpy as np
import pandas as pd
from scipy.stats import norm

from . import _constants as C
from .report import CalibrationRecord, CIResult, CalibratedResult

logger = logging.getLogger("openreltime")

__all__ = ["confidence_interval"]

#: relative tolerance used to decide that a node age sits on a hard bound
_PIN_TOL = 1e-9


@dataclass
class BranchVariance:
    """Sampling variance ``vS(b)`` keyed by the child node id of each edge."""

    values: dict[int, float]


def _central_derivative_column(
    func, base_vector: np.ndarray, index: int
) -> np.ndarray:
    """Column of partial derivatives d(func)/d(b_index) by central differences."""
    h = max(1e-8, 1e-5 * abs(base_vector[index]))
    up = base_vector.copy()
    up[index] += h
    dn = base_vector.copy()
    dn[index] -= h
    return (func(up) - func(dn)) / (2 * h)


def _pipeline_settings(calibrated: CalibratedResult) -> dict[str, Any]:
    """Estimator settings the differentiation pipeline has to reproduce."""
    params = dict(calibrated.params or {})
    mean = params.get("mean", "geometric")
    if mean not in ("geometric", "arithmetic"):
        raise ValueError(
            f"calibrated.params['mean']={mean!r} is not a known averaging "
            "convention; confidence_interval can only differentiate the point "
            "estimate if it reproduces it ('geometric' or 'arithmetic')"
        )
    method = params.get("method", "bounds")
    if method not in ("bounds", "effective"):
        raise ValueError(
            f"confidence_interval supports calibrated results computed with "
            f"method='bounds' or method='effective', got {method!r}"
        )
    return {
        "mean": mean,
        "method": method,
        # a missing key means the default guard, an explicit None means "off"
        "rate_ratio_threshold": params.get(
            "rate_ratio_threshold", C.RATE_RATIO_THRESHOLD
        ),
    }


def _pipeline_calibrations(
    calibrated: CalibratedResult, settings: dict[str, Any]
) -> list[Any]:
    """The calibration rows the point estimate was actually computed with.

    ``method="effective"`` results store the effective bounds in
    ``CalibrationRecord.used_min_bound/used_max_bound`` (msz236 uses effective
    bounds *as* calibration constraints), so the derivatives follow the
    effective-bounds path instead of the raw - possibly density-only, hence
    unbounded - rows the user supplied.
    """
    from .calibrate import Calibration

    records = calibrated.calibrations
    if settings["method"] == "effective":
        rows = [
            Calibration(
                node_id=rec.node_id,
                min_bound=rec.used_min_bound,
                max_bound=rec.used_max_bound,
            )
            for rec in records
            if rec.used_min_bound is not None or rec.used_max_bound is not None
        ]
        if not rows:
            raise ValueError(
                "the calibrated result carries no effective bounds, so the CI "
                "pipeline cannot reproduce its point estimate; re-run "
                "calibrate(method='effective') with n_effective >= 2, or use "
                "the effective_bounds table directly"
            )
        return rows
    return [
        Calibration(node_id=rec.node_id, min_bound=rec.min_bound, max_bound=rec.max_bound)
        for rec in records
    ]


def _make_pipeline(
    calibrated: CalibratedResult, kind: str
) -> Callable[[np.ndarray], np.ndarray]:
    """Build a function branch-vector -> node vector for one estimator setting.

    ``kind="rates"`` returns the RRF relative lineage rates ``r_j`` of eqs
    (1)-(4).  They are functions of the branch lengths only - no calibration
    enters - which is exactly what eqs (12)-(13) differentiate.
    ``kind="times"`` returns the absolute node times along the full estimation
    path (RRF times -> global factor -> hard-bound adjustment), which is what
    eq (7) differentiates.  ``kind="relative"`` returns ``f * t_rel``, i.e. the
    same path with the hard-bound adjustment skipped; it is only used for the
    fallback sensitivity of nodes pinned at a one-sided bound.

    The trial tree is a copy, so the caller's tree is never perturbed, and the
    returned vectors are indexed by the *trial* tree's nodes after an explicit
    node-id alignment check.
    """
    from .calibrate import calibrate as _calibrate
    from .times import rrf_times as _rrf_times

    if kind not in ("rates", "times", "relative"):
        raise ValueError(f"unknown pipeline kind {kind!r}")
    settings = _pipeline_settings(calibrated)
    calibs = _pipeline_calibrations(calibrated, settings) if kind == "times" else []
    f0 = float(calibrated.time_factor)
    tree = calibrated.tree
    base_edge_ids = [n.node_id for n in tree.walk() if not n.is_root()]
    base_internal_ids = [n.node_id for n in tree.walk() if not n.is_tip()]

    def evaluate(b_vector: np.ndarray) -> np.ndarray:
        trial = tree.copy()
        trial_edges = [n for n in trial.walk() if not n.is_root()]
        if len(trial_edges) != len(b_vector):
            raise ValueError(
                f"branch vector of length {len(b_vector)} does not match the "
                f"{len(trial_edges)} branches of the analysed tree"
            )
        for i, n in enumerate(trial_edges):
            n.blen = float(b_vector[i])
        trial.assign_ids()
        trial_edge_ids = [n.node_id for n in trial.walk() if not n.is_root()]
        trial_internal_ids = [n.node_id for n in trial.walk() if not n.is_tip()]
        if trial_edge_ids != base_edge_ids or trial_internal_ids != base_internal_ids:
            raise ValueError(
                "node-id space of the perturbed trial tree diverged from the "
                "analysed tree; the CI derivatives would be indexed against the "
                f"wrong nodes ({trial_internal_ids[:5]}... != "
                f"{base_internal_ids[:5]}...)"
            )
        res = _rrf_times(
            trial,
            mean=settings["mean"],
            rate_ratio_threshold=settings["rate_ratio_threshold"],
        )
        if kind == "rates":
            return np.array(
                [res.rates[n.node_id] for n in trial.walk() if not n.is_root()],
                dtype=float,
            )
        if kind == "relative":
            return np.array(
                [f0 * res.times[n.node_id] for n in trial.walk() if not n.is_tip()],
                dtype=float,
            )
        if calibs:
            cal = _calibrate(res, calibs, method="bounds")
            times_map = cal.times
        else:
            times_map = res.times
        return np.array(
            [times_map[n.node_id] for n in trial.walk() if not n.is_tip()],
            dtype=float,
        )

    return evaluate


def _calibration_uncertainty(
    rec: CalibrationRecord,
) -> tuple[Optional[float], Optional[float], str]:
    """The interval a calibration itself asserts for the node age.

    Preference order: effective bounds (msz236 2.5/97.5 percentiles of the
    replicated inferences) -> the 2.5/97.5 quantiles of a supplied density ->
    the user-supplied hard bounds.  ``None`` on a side means the calibration is
    one-sided there.
    """
    hard = (rec.min_bound, rec.max_bound)
    used = (rec.used_min_bound, rec.used_max_bound)
    if used != (None, None) and used != hard:
        lo, hi = used
        return lo, hi, "effective bounds (msz236 2.5/97.5 percentiles)"
    if rec.density and rec.density_params:
        from .calibrate import density_quantiles

        try:
            lo, hi = density_quantiles(rec.density, rec.density_params)
        except ValueError:  # pragma: no cover - guarded at parse time
            lo, hi = float("nan"), float("nan")
        if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
            if rec.min_bound is not None:
                lo = max(lo, rec.min_bound)
            if rec.max_bound is not None:
                hi = min(hi, rec.max_bound)
            if hi > lo:
                return lo, hi, f"{rec.density} density quantiles"
    if rec.min_bound is not None and rec.max_bound is not None:
        return rec.min_bound, rec.max_bound, "user-supplied min-max interval"
    return rec.min_bound, rec.max_bound, "one-sided hard bound"


def _free_se_array(
    calibrated: CalibratedResult,
    blens: np.ndarray,
    v_b: np.ndarray,
    n_internal: int,
) -> np.ndarray:
    """sqrt of eq (7) evaluated on the *unconstrained* path ``t = f * t_rel``.

    Only needed when a node is pinned at a one-sided hard bound, where the
    calibrated pipeline returns an exactly zero derivative.
    """
    pipeline = _make_pipeline(calibrated, "relative")
    derivs = np.zeros((n_internal, blens.size))
    for i in range(blens.size):
        if v_b[i] == 0:
            continue
        derivs[:, i] = _central_derivative_column(pipeline, blens, i)
    return np.sqrt(np.sum(derivs**2 * v_b[None, :], axis=1))


def confidence_interval(
    calibrated: CalibratedResult,
    *,
    branch_var: Optional[dict[int, float]] = None,
    n_sites: Optional[int] = None,
    seq_length: Optional[int] = None,
    level: float = 0.95,
    outgroup_check: Optional[str] = None,
    **params: Any,
) -> CIResult:
    """Analytical confidence interval for calibrated node times (msz236).

    Parameters
    ----------
    calibrated:
        Result of :func:`openreltime.calibrate`.
    branch_var:
        User-provided sampling variance ``vS(b)`` keyed by child node id
        (highest priority source).
    n_sites / seq_length:
        Poisson approximation ``vS(b) = b / L`` when no explicit variances
        are given (documented approximation).
        When neither is supplied, ``vS = 0`` and the interval captures only
        the rate-heterogeneity component (the msz236 simulation protocol).
    level:
        Confidence level (default 0.95).
    outgroup_check:
        Provenance of the rooting check performed when the tree was read
        (``"error"`` / ``"warn"``); echoed into ``CIResult.params`` together
        with the value recorded by the calibration step, so a CI table can be
        traced back to the exact input handling that produced it.
    **params:
        Any further provenance keys are recorded verbatim in
        ``CIResult.params``.

    Returns
    -------
    CIResult
        Per-node table with standard error, truncated bounds and the
        ``se_reliable``/``notes`` flags.
    """
    tree = calibrated.tree
    internal_nodes = [n for n in tree.walk() if not n.is_tip()]
    edge_nodes = [n for n in tree.walk() if not n.is_root()]
    blens = np.array([float(n.blen or 0.0) for n in edge_nodes])
    n_branches = blens.size  # N of msz236 eqs (7)-(14)
    if outgroup_check is None:
        outgroup_check = (calibrated.params or {}).get("outgroup_check")

    # -- vS(b) ------------------------------------------------------------------
    if branch_var:
        v_s = np.array([float(branch_var.get(n.node_id, 0.0)) for n in edge_nodes])
        v_s_source = "user-provided branch variance"
    elif seq_length:
        v_s = blens / float(seq_length)
        v_s_source = f"Poisson approximation b/L (L={seq_length})"
    elif n_sites:
        v_s = blens / float(n_sites)
        v_s_source = f"Poisson approximation b/L (L={n_sites})"
    else:
        v_s = np.zeros_like(blens)
        v_s_source = "zero (rate-heterogeneity variance only; msz236 protocol)"

    settings = _pipeline_settings(calibrated)

    # -- Vobs(R), eq (9), on the RRF relative rates r_j of eqs (1)-(4) ----------
    # calibrated.rates are the RRF rates; the implied absolute rates of the
    # calibrated tree live in calibrated.implied_rates and are a different
    # physical quantity.  Lineages without a finite rate (e.g. a zero-length
    # calibrated duration) are dropped from the average instead of entering it
    # as a spurious 0.0.
    rates = np.array(
        [calibrated.rates.get(n.node_id, float("nan")) for n in edge_nodes],
        dtype=float,
    )
    usable = np.isfinite(rates)
    n_lineages = int(usable.sum())
    n_dropped = int(n_branches - n_lineages)
    if n_dropped:
        logger.warning(
            "%d of %d lineages have no finite RRF rate and are dropped from "
            "Vobs(R) and SV(R)",
            n_dropped,
            n_branches,
        )
    if n_lineages < 2:
        raise ValueError(
            f"only {n_lineages} lineage(s) carry a finite RRF rate: Vobs(R) "
            "(msz236 eq 9) cannot be estimated. Check the branch lengths and the "
            "rate-ratio guard of the underlying rrf_times run."
        )
    r_used = rates[usable]
    observed_var = float(np.mean((r_used - r_used.mean()) ** 2))  # eq (9)

    # -- SV(R) (eqs 11-12) and RV(R) (eqs 10, 13) ------------------------------
    # sv(r_j) = sum_i (dr_j/db_i)^2 vS(b_i)  (eq 12);
    # SV(R)   = sum_j sv(r_j)                (eq 11), computed here as the
    #           branch-wise accumulation below;
    # RV(R)   = Vobs(R) - SV(R)/N             (eqs 10, 13 on the per-lineage
    #           scale of eq 9 - eq 13 as printed subtracts the *sum* SV(R) from
    #           the *average* Vobs(R), which clamps RV(R) to 0 for every
    #           non-trivial tree).
    sv_total = 0.0
    if v_s.any():
        rate_pipeline = _make_pipeline(calibrated, "rates")
        for i in range(n_branches):
            if v_s[i] == 0:
                continue
            column = _central_derivative_column(rate_pipeline, blens, i)[usable]
            sv_total += float(np.sum(column**2)) * v_s[i]
    sv_per_lineage = sv_total / n_lineages
    rv = max(observed_var - sv_per_lineage, 0.0)
    if sv_per_lineage > observed_var:
        logger.warning(
            "sampling variance exceeds observed rate variance "
            "(SV(R)/N=%.3g > Vobs(R)=%.3g, i.e. SV(R)=%.3g over %d lineages); "
            "RV(R) clamped to 0 and the interval rests on vS(b) alone",
            sv_per_lineage,
            observed_var,
            sv_total,
            n_lineages,
        )

    # -- vR(b) (eq 14) and total branch variance (eq 8) --------------------------
    denom = float(np.sum(blens**2))
    v_r = (blens**2 / denom) * rv if denom > 0 else np.zeros_like(blens)
    v_b = v_s + v_r

    # -- node-time variance (eq 7 via whole-pipeline derivatives) ----------------
    time_pipeline = _make_pipeline(calibrated, "times")
    derivs = np.zeros((len(internal_nodes), n_branches))
    for i in range(n_branches):
        if v_b[i] == 0:
            continue
        derivs[:, i] = _central_derivative_column(time_pipeline, blens, i)
    v_t = np.sum(derivs**2 * v_b[None, :], axis=1)

    # -- interval with hard-bound truncation --------------------------------------
    t_abs = np.array([calibrated.times[n.node_id] for n in internal_nodes], dtype=float)
    se = np.sqrt(v_t)
    z = float(norm.ppf(0.5 + level / 2))
    lower = t_abs - z * se
    upper = t_abs + z * se
    reliable = np.ones(len(internal_nodes), dtype=bool)
    notes: list[str] = [""] * len(internal_nodes)
    pos = {n.node_id: i for i, n in enumerate(internal_nodes)}
    free_se: Optional[np.ndarray] = None

    for rec in calibrated.calibrations:
        idx = pos.get(rec.node_id)
        if idx is None:
            continue
        t = float(t_abs[idx])
        scale = max(1.0, abs(t))
        # msz236: truncate at the imposed hard calibration constraints
        truncated = False
        if rec.min_bound is not None and lower[idx] < rec.min_bound:
            lower[idx] = rec.min_bound
            truncated = True
        if rec.max_bound is not None and upper[idx] > rec.max_bound:
            upper[idx] = rec.max_bound
            truncated = True
        if lower[idx] > t or upper[idx] < t:
            # the solver could not reach this constraint (it records a warning);
            # an inverted interval must not be printed, and the node is not
            # qualified by eq (7) either
            lower[idx] = min(lower[idx], t)
            upper[idx] = max(upper[idx], t)
            notes[idx] = "calibration constraint not satisfied by the solver; "
            reliable[idx] = False
        elif truncated:
            notes[idx] = "truncated at a hard calibration bound; "

        at_min = rec.min_bound is not None and abs(t - rec.min_bound) <= _PIN_TOL * scale
        at_max = rec.max_bound is not None and abs(t - rec.max_bound) <= _PIN_TOL * scale
        if se[idx] > _PIN_TOL * scale:
            continue
        # The age of a calibrated node is a fixed point of the calibration: either
        # the bound itself (an attractor of the violation loop, so t(b+h) = t(b-h))
        # or the midpoint of the feasible interval (t = (min+max)/2 whatever the
        # branch lengths).  eq (7) therefore returns exactly zero for it, which
        # must not be printed as "zero uncertainty in age": the node inherits
        # the uncertainty its own calibration asserts instead.
        lo_cal, hi_cal, source = _calibration_uncertainty(rec)
        how = (
            "pinned at its minimum hard bound"
            if at_min
            else "pinned at its maximum hard bound"
            if at_max
            else "set by the calibration interval"
        )
        cand_lo = cand_hi = None
        if lo_cal is not None and hi_cal is not None:
            cand_lo, cand_hi = min(lo_cal, t), max(hi_cal, t)
            # msz236: the reported interval never leaves the imposed constraints
            if rec.min_bound is not None:
                cand_lo = max(cand_lo, rec.min_bound)
            if rec.max_bound is not None:
                cand_hi = min(cand_hi, rec.max_bound)
        if cand_lo is not None and cand_hi - cand_lo > _PIN_TOL * scale:
            lower[idx], upper[idx] = cand_lo, cand_hi
        else:
            # the calibration asserts no width here (a one-sided bound, or an
            # effective interval collapsed by the hard bound): the open side
            # comes from the sensitivity of the unconstrained path t = f * t_rel
            if free_se is None:
                free_se = _free_se_array(calibrated, blens, v_b, len(internal_nodes))
            delta = z * float(free_se[idx])
            if delta <= 0:
                notes[idx] += (
                    f"age {how} at a one-sided hard bound and no branch variance "
                    "is available; interval degenerate"
                )
                reliable[idx] = False
                continue
            if at_max:
                lower[idx], upper[idx] = t - delta, t
            elif at_min:
                lower[idx], upper[idx] = t, t + delta
            else:
                lower[idx], upper[idx] = t - delta, t + delta
            source += ", open side from the unconstrained path f*t_rel"
        se[idx] = float(upper[idx] - lower[idx]) / (2 * z) if z > 0 else se[idx]
        notes[idx] += (
            f"age {how}: interval = the calibration's own {source}, not the "
            "eq-(7) delta method"
        )
        reliable[idx] = False

    for i, node in enumerate(internal_nodes):
        if lower[i] < 0.0:
            # a divergence time cannot be negative; the symmetric delta-method
            # interval is truncated at the origin for very young nodes
            lower[i] = 0.0
            notes[i] += "lower bound clipped at 0 (ages are non-negative); "
        if reliable[i] and se[i] <= _PIN_TOL * max(1.0, abs(t_abs[i])):
            notes[i] += (
                "zero delta-method variance (no branch variance: vS(b)=0 and "
                "RV(R)=0)"
            )
            reliable[i] = False
    if not reliable.any():
        logger.warning(
            "no node time is qualified by the eq-(7) branch-variance "
            "decomposition in this analysis; see the notes column"
        )

    table = pd.DataFrame(
        {
            "node_id": [n.node_id for n in internal_nodes],
            "label": [n.label or "-" for n in internal_nodes],
            "time": t_abs,
            "se": se,
            "lower": lower,
            "upper": upper,
            "width": upper - lower,
            "se_reliable": reliable,
            "notes": [n.rstrip("; ") for n in notes],
        }
    )
    method_notes = (
        "delta method (msz236 eqs 7-14); vS source: "
        f"{v_s_source}; Vobs(R)={observed_var:.4g}; SV(R)/N={sv_per_lineage:.4g}; "
        f"RV(R)={rv:.4g}; pipeline: method={settings['method']}, "
        f"mean={settings['mean']}"
    )
    merged = {
        "level": level,
        "v_s_source": v_s_source,
        "outgroup_check": outgroup_check,
        "Vobs_R": observed_var,
        "SV_R": sv_total,
        "SV_per_lineage_R": sv_per_lineage,
        "RV_R": rv,
        "n_branches": int(n_branches),
        "n_lineages": n_lineages,
        "n_lineages_dropped": n_dropped,
        "rate_definition": "RRF relative lineage rates (msz236 eqs 1-4)",
        "pipeline_method": settings["method"],
        "pipeline_mean": settings["mean"],
        "n_nodes_se_not_reliable": int((~reliable).sum()),
        **params,
    }
    return CIResult(table=table, level=level, method=method_notes, params=merged)
