"""Calibration of relative times and absolute-time conversion.

Implements the calibration mechanism of MEGA RelTime as described in

    Tao Q, Tamura K, Mello B, Kumar S.  *Reliable confidence intervals for
    RelTime estimates of evolutionary divergence times.*  Mol Biol Evol
    37:280-290 (2020), doi:10.1093/molbev/msz236

* ``bounds`` method: the global time factor ``f`` is searched as the midpoint
  of the intersection of the per-calibration feasible intervals
  ``[min/t, max/t]``; a one-sided bound set takes its finite edge instead, which
  keeps the result invariant to the unit of the supplied ages.  Residual
  violations trigger proportional rate scaling of the offending lineage
  (hard-bound semantics of Tamura et al. 2013), propagated to descendants,
  iterated up to ``CALIBRATION_MAX_ITER`` times.
* ``effective`` method (msz236): for calibrations carrying a probability
  density, two dates are sampled from the density per replicate and used as
  a (min, max) pair; the analysis is repeated ``n_effective`` times and the
  2.5/97.5 percentiles of the resulting node-age distributions become
  *effective bounds*, which are then used for the final point estimate.

The fine iteration order of MEGA's violation propagation is a heuristic not
fully specified by the paper; this implementation follows the msz236 prose
(and the ``TimeFactor``/``MinTimeFactor``/``MaxTimeFactor`` design of MEGA's
``mreltimecomputer.pas``, GPL-3) and is validated by tests against
hand-worked cases (Gate G7, docs/validation.md).

Input validation is strict on purpose: a calibration row that carries neither
a finite ``min_bound`` nor a finite ``max_bound`` contributes no constraint at
all, so the feasible-factor search would return ``f = inf`` and every absolute
age would silently become ``inf``.  Such rows are rejected at
construction time, and ``method="bounds"`` additionally rejects density-only
rows, which only the effective-bounds procedure can use.  ``density_params`` is
checked against a per-family key whitelist (:data:`DENSITY_PARAM_KEYS`) so a
mis-spelled prior cannot be silently sampled from the family's defaults, and an
unsatisfiable nested bound (a descendant older than its ancestor) raises instead
of returning a result that violates the caller's own bounds.  That guarantee is
checked twice: once on the worst violator inside the repair loop and once on the
converged vector as a whole, so a solution that ends non-monotonic or with a
negative age is rejected rather than written out with a warning.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Optional

import numpy as np
import pandas as pd
from scipy import stats

from . import _constants as C
from .report import CalibrationRecord, CalibratedResult
from .taxonomy import resolve_taxa_tokens, warn_if_not_monophyletic
from .tree import PhyloNode

logger = logging.getLogger("openreltime")

__all__ = [
    "Calibration",
    "CalibrationCancelled",
    "DENSITY_PARAM_KEYS",
    "DENSITY_TYPES",
    "density_quantiles",
    "parse_calibrations",
    "calibrate",
]

DENSITY_TYPES = ("uniform", "exponential", "normal", "lognormal")

#: The *complete* parameter key set of every supported density family.
#:
#: ``_sample_density_pair`` and :func:`density_quantiles` used to read each key
#: through ``params.get(key, default)``, so ``density='normal',
#: density_params='mu=5;sigma=1'`` was silently sampled as ``N(0, 1)``: the
#: intended prior was discarded and a different one substituted without a word.
#: Keys are therefore validated against this whitelist, which is also what makes
#: the documented "parameters are parsed strictly" claim in ``docs/parameters.md``
#: true: an unknown key (a typo, or a different family's vocabulary) and a
#: missing required key both raise instead of being ignored or defaulted.
DENSITY_PARAM_KEYS: dict[str, tuple[str, ...]] = {
    "uniform": ("min", "max"),
    "exponential": ("offset", "mean"),
    "normal": ("mean", "sd"),
    "lognormal": ("offset", "meanlog", "sdlog"),
}


class CalibrationCancelled(RuntimeError):
    """Raised when a caller-requested cancellation arrives mid-calibration."""


#: Relative slack of the post-solve time-order check in
#: :meth:`_CalibratedTree._check_time_order`.  ``_scale_subtree`` rewrites every
#: duration of a lineage in floating point, so a legitimate solution can put a
#: descendant a few ULP above its ancestor; only a real inversion is rejected.
#: Same order as the ``_PIN_TOL`` the CI layer uses for "sits on the bound".
_MONOTONIC_TOL = 1e-9


def _fmt_bound(value: Optional[float]) -> str:
    """Render one calibration bound for an error message (``None`` = open)."""
    return "none" if value is None else f"{value:.6g}"


def _validate_density_params(
    density: str, params: dict[str, float], *, context: str
) -> None:
    """Raise unless *params* holds exactly the keys *density* declares.

    ``context`` prefixes the message with whatever identifies the offending
    row (a ``node_id``/``taxon_set`` description, or a file and line number).
    """
    try:
        required = DENSITY_PARAM_KEYS[density]
    except KeyError:  # pragma: no cover - density family checked at construction
        raise ValueError(
            f"{context}: unsupported density {density!r}; "
            f"choose from {DENSITY_TYPES}"
        ) from None
    unknown = sorted(set(params) - set(required))
    if unknown:
        raise ValueError(
            f"{context}: density {density!r} accepts only the keys "
            f"{list(required)}, but {unknown} were given. Those keys are "
            "ignored by no sampler: they are rejected, because silently "
            "defaulting them would substitute a different prior for the one "
            "you asked for."
        )
    missing = [key for key in required if key not in params]
    if missing:
        raise ValueError(
            f"{context}: density {density!r} needs the keys "
            f"{list(required)}, but {missing} are missing. Every parameter of "
            "the family must be stated explicitly; none is defaulted."
        )


@dataclass(frozen=True)
class Calibration:
    """A calibration constraint.

    Either ``node_id`` or ``taxon_set`` identifies the target node.  A row must
    carry at least one *finite* bound (``min_bound`` / ``max_bound``) or a
    ``density``: a row with neither constrains nothing, and under
    ``method="bounds"`` densities are not used at all, so such a row is rejected
    here rather than being allowed to drive the global time factor to ``inf``.
    ``density_params`` keys: uniform ``min,max``; exponential ``offset,mean``;
    normal ``mean,sd``; lognormal ``offset,meanlog,sdlog``.  The key set is
    validated against :data:`DENSITY_PARAM_KEYS`: an unknown key or a missing
    required key raises, so a mis-spelled prior can never be silently sampled
    from the family's defaults.
    """

    taxon_set: Optional[frozenset[str]] = None
    node_id: Optional[int] = None
    min_bound: Optional[float] = None
    max_bound: Optional[float] = None
    density: Optional[str] = None
    density_params: dict[str, float] = field(default_factory=dict)

    @property
    def has_bound(self) -> bool:
        """True when the row carries at least one bound."""
        return self.min_bound is not None or self.max_bound is not None

    @property
    def is_density_only(self) -> bool:
        """True when the row carries a density but no bound (effective only)."""
        return self.density is not None and not self.has_bound

    def describe(self) -> str:
        """Human-readable identity of the target, for error messages."""
        if self.node_id is not None:
            return f"node_id={self.node_id}"
        if self.taxon_set:
            return "taxon_set=" + "|".join(sorted(self.taxon_set))
        return "<calibration without node_id or taxon_set>"

    def __post_init__(self) -> None:
        for name, value in (("min_bound", self.min_bound), ("max_bound", self.max_bound)):
            if value is None:
                continue
            if not math.isfinite(value):
                raise ValueError(
                    f"{self.describe()}: {name}={value!r} is not a finite number; "
                    "calibration bounds must be finite ages"
                )
            if value < 0:
                raise ValueError(
                    f"{self.describe()}: {name}={value!r} is negative; node ages "
                    "and therefore calibration bounds are non-negative"
                )
        if self.max_bound is not None and self.max_bound <= 0:
            raise ValueError(
                f"{self.describe()}: max_bound must be > 0 (a zero maximum forces "
                "every absolute age to zero)"
            )
        if (
            self.min_bound is not None
            and self.max_bound is not None
            and self.min_bound > self.max_bound
        ):
            raise ValueError(
                f"{self.describe()}: min_bound={self.min_bound} exceeds "
                f"max_bound={self.max_bound}; the calibration interval is empty"
            )
        if self.density is not None and self.density not in DENSITY_TYPES:
            raise ValueError(
                f"{self.describe()}: unsupported density {self.density!r}; "
                f"choose from {DENSITY_TYPES}"
            )
        if self.density is not None:
            _validate_density_params(
                self.density, self.density_params, context=self.describe()
            )
        if not self.has_bound and self.density is None:
            raise ValueError(
                f"{self.describe()}: a calibration needs at least one finite bound "
                "(min_bound and/or max_bound) or a probability density, but "
                "neither was given.  Under method='bounds' a density does not "
                "count as a bound: add min_bound/max_bound, or run "
                "calibrate(..., method='effective') for density-only rows."
            )


def parse_calibrations(path: str | Path) -> list[Calibration]:
    """Read a ``calibrations.tsv`` file.

    Columns: ``node_id`` (optional), ``taxon_set`` (``|``-separated, optional;
    entries may be tip labels *or* group names resolved against the tree via
    ``--taxon-table``/auto-detected label prefixes),
    ``min_bound``, ``max_bound`` (``.`` = missing), ``density`` (optional),
    ``density_params`` (semicolon-separated ``key=value`` pairs, optional).

    Only explicit ``key=value`` parsing is performed; expression evaluation
    (``eval``) is deliberately never used on this input.  ``density_params`` is
    read *strictly*: every key must belong to the declared family's key set
    (:data:`DENSITY_PARAM_KEYS`), all of the family's keys must be present, and a
    non-empty ``density_params`` on a row without a ``density`` raises rather
    than being dropped.
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"min_bound", "max_bound"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"calibration file missing columns: {sorted(missing)!r}")

    def _num(value: str) -> Optional[float]:
        value = value.strip()
        if value in ("", ".", "NA", "-"):
            return None
        try:
            return float(value)
        except ValueError as exc:
            raise ValueError(f"cannot parse numeric bound {value!r} in {path}") from exc

    calibrations: list[Calibration] = []
    for line, (_, row) in enumerate(df.iterrows(), start=2):
        node_id: Optional[int] = None
        taxon_set: Optional[frozenset[str]] = None
        if str(row.get("node_id", "")).strip() not in ("", ".", "NA", "-"):
            node_id = int(str(row["node_id"]).strip())
        taxa = str(row.get("taxon_set", "")).strip()
        if taxa not in ("", ".", "NA", "-"):
            taxon_set = frozenset(t for t in taxa.split("|") if t)
        if node_id is None and taxon_set is None:
            raise ValueError(
                f"{path} line {line}: calibration row without node_id or taxon_set: "
                f"{row.to_dict()}"
            )
        density = str(row.get("density", "")).strip().lower()
        density = density if density not in ("", ".", "NA", "-") else None
        if density is not None and density not in DENSITY_TYPES:
            raise ValueError(
                f"{path} line {line}: unsupported density {density!r}; "
                f"choose from {DENSITY_TYPES}"
            )
        dparams: dict[str, float] = {}
        raw_params = str(row.get("density_params", "")).strip()
        if raw_params not in ("", ".", "NA", "-"):
            if density is None:
                raise ValueError(
                    f"{path} line {line}: density_params={raw_params!r} given "
                    "without a density column; there is no family to read the "
                    "keys against, so they would be silently ignored -- add "
                    f"density (one of {DENSITY_TYPES}) or drop density_params"
                )
            for chunk in raw_params.split(";"):
                chunk = chunk.strip()
                if not chunk:
                    continue
                if "=" not in chunk:
                    raise ValueError(
                        f"{path} line {line}: malformed density_params chunk "
                        f"{chunk!r}; expected key=value"
                    )
                key, _, value = chunk.partition("=")
                dparams[key.strip()] = float(value.strip())
        try:
            calibration = Calibration(
                taxon_set=taxon_set,
                node_id=node_id,
                min_bound=_num(str(row["min_bound"])),
                max_bound=_num(str(row["max_bound"])),
                density=density,
                density_params=dparams,
            )
        except ValueError as exc:
            # an unconstrained row must be reported where it is read, not
            # silently turned into f = inf later on
            raise ValueError(f"{path} line {line}: {exc}") from exc
        calibrations.append(calibration)
    if not calibrations:
        raise ValueError(f"no calibrations found in {path}")
    return calibrations


def _sample_density_pair(cal: Calibration, rng: np.random.Generator) -> tuple[float, float]:
    """Sample two ordered dates from a calibration density (msz236).

    The keys are read directly, never through ``dict.get`` with a fallback: the
    family's parameter set was already validated when *cal* was constructed, so
    a default here would be a silent mis-specification rather than a convenience.
    """
    p = cal.density_params
    assert cal.density is not None
    if cal.density == "uniform":
        a, b = p["min"], p["max"]
        x = rng.uniform(a, b, size=2)
    elif cal.density == "exponential":
        offset, mean = p["offset"], p["mean"]
        x = offset + rng.exponential(mean, size=2)
    elif cal.density == "normal":
        mu, sd = p["mean"], p["sd"]
        x = rng.normal(mu, sd, size=2)
    elif cal.density == "lognormal":
        offset = p["offset"]
        mu, sd = p["meanlog"], p["sdlog"]
        x = offset + rng.lognormal(mu, sd, size=2)
    else:  # pragma: no cover - validated at parse time
        raise ValueError(cal.density)
    return float(min(x)), float(max(x))


def density_quantiles(
    density: str,
    params: dict[str, float],
    probs: tuple[float, float] = (0.025, 0.975),
) -> tuple[float, float]:
    """Two quantiles of a supported calibration density.

    msz236 summarises a calibration density by its 2.5/97.5 percentiles (the
    effective bounds); :mod:`openreltime.ci` uses the same pair as the
    *calibration's own* uncertainty for nodes whose age is pinned at a hard
    bound, where the delta method of eq (7) returns exactly zero.

    ``params`` is validated against :data:`DENSITY_PARAM_KEYS` before it is
    used, so a partial or mis-spelled parameter set raises instead of being
    filled in from the family's defaults.
    """
    _validate_density_params(
        density, params, context=f"density_quantiles({density!r})"
    )
    if density == "uniform":
        a, b = params["min"], params["max"]
        dist = stats.uniform(loc=a, scale=max(b - a, 0.0))
    elif density == "exponential":
        dist = stats.expon(loc=params["offset"], scale=params["mean"])
    elif density == "normal":
        dist = stats.norm(loc=params["mean"], scale=params["sd"])
    elif density == "lognormal":
        dist = stats.lognorm(
            s=params["sdlog"],
            scale=math.exp(params["meanlog"]),
            loc=params["offset"],
        )
    else:
        raise ValueError(f"unsupported density {density!r}; choose from {DENSITY_TYPES}")
    low, high = (float(q) for q in dist.ppf(probs))
    return low, high


class _CalibratedTree:
    """Absolute-time computation for a fixed set of (min, max) bounds."""

    def __init__(
        self,
        tree: PhyloNode,
        rel_result,
        bounds: list[tuple[PhyloNode, Optional[float], Optional[float]]],
    ) -> None:
        self.tree = tree
        self.rel_result = rel_result
        self.bounds = bounds
        self.rel_times = np.array(
            [rel_result.times[n.node_id] for n in tree.walk() if not n.is_tip()]
        )
        self.internal_nodes = [n for n in tree.walk() if not n.is_tip()]
        self.pos = {n.node_id: i for i, n in enumerate(self.internal_nodes)}
        self.root_pos = self.pos[tree.node_id]

    # -- core ---------------------------------------------------------------

    def solve(self) -> tuple[float, np.ndarray, list[str]]:
        """Return (f, absolute node times, log messages) for the fixed bounds.

        Raises :class:`ValueError` if the converged ages are not nested in the
        tree (see :meth:`_check_time_order`): a solution the violation loop
        reports as satisfied is still nonsense if it puts a descendant older
        than its ancestor or below the present.
        """
        times = self.rel_times.copy()
        messages: list[str] = []
        f = self._feasible_factor(times)
        times = times * f
        times, iter_log = self._enforce_bounds(times)
        messages.extend(iter_log)
        self._check_time_order(times)
        if self.root_pos < times.size and times[self.root_pos] > 0:
            f = float(times[self.root_pos] / max(self.rel_times[self.root_pos], 1e-300))
        return f, times, messages

    def _feasible_factor(self, times: np.ndarray) -> float:
        """Midpoint of ``[max(min/t), min(max/t)]`` (msz236, Tamura et al. 2013).

        Every row must contribute a finite bound: a row with neither ``min`` nor
        ``max`` leaves the interval unbounded on both sides, which would return
        ``inf`` and turn every absolute age into ``inf`` without a word of
        complaint.

        A one-sided bound set takes its finite edge: with only minimum bounds
        ``f = max(min/t)``, with only maximum bounds ``f = min(max/t)``, so the
        binding node lands exactly on that bound.  Both edges live on the same
        relative scale as the ``min``/``max`` ages the caller supplied, so the
        solution is invariant to the unit of those ages.  (An earlier release
        floored either edge at ``1.0``, a dimensionless constant of the
        arbitrary RRF relative scale: it left a minimum bound that asked for a
        *younger* tree inert and inflated a maximum bound above it, and it made
        the same physical constraint behave differently in Ma and in Ga.)
        """
        lo, hi = -np.inf, np.inf
        unbounded = [
            node.node_id for node, mn, mx in self.bounds if mn is None and mx is None
        ]
        if unbounded:
            raise ValueError(
                "no finite calibration bound for node(s) "
                f"{sorted(unbounded)}: the 'bounds' solver cannot derive a global "
                "time factor from rows that carry neither min_bound nor max_bound "
                "(a density alone is only usable with method='effective').  Add a "
                "min_bound and/or max_bound to those rows, or drop them."
            )
        for node, mn, mx in self.bounds:
            t = max(times[self.pos[node.node_id]], 1e-300)
            if mn is not None:
                lo = max(lo, mn / t)
            if mx is not None:
                hi = min(hi, mx / t)
        if lo > hi:
            # infeasible: start from the midpoint of the overlap gap so the
            # violation loop can migrate the solution toward feasibility
            f = (lo + hi) / 2 if math.isfinite(lo) and math.isfinite(hi) else max(lo, hi)
        elif not math.isfinite(lo):
            f = hi  # max-only bound set: the finite edge itself
        elif not math.isfinite(hi):
            f = lo  # min-only bound set: the finite edge itself
        else:
            f = (lo + hi) / 2
        if not math.isfinite(f) or f <= 0:
            raise ValueError(
                f"the calibration set admits no usable time factor (f={f!r}): "
                "check that every relative node time of the calibrated nodes is "
                "positive and that the bounds are mutually compatible (a bound "
                "set whose only finite edge is 0 fixes no time scale) -> "
                + "; ".join(
                    f"node {node.node_id} [{mn}, {mx}]"
                    for node, mn, mx in self.bounds
                )
            )
        return float(f)

    def _enforce_bounds(self, times: np.ndarray) -> tuple[np.ndarray, list[str]]:
        """Hard-bound violation adjustment (msz236 / Tamura et al. 2013).

        A violation the geometry of the tree cannot absorb - a descendant bound
        whose ancestor is already younger than the bound asks it to be - raises
        :class:`ValueError` rather than returning: the remaining rows would be
        abandoned mid-loop and the returned ages would violate the very
        ``min_bound``/``max_bound`` the caller supplied.

        The nested check below covers only the *worst violator of the current
        round*, so it cannot see an inversion the loop has just created
        elsewhere in the tree; :meth:`_check_time_order` re-checks the whole
        converged vector, which is what makes the guarantee total rather than
        per-round.
        """
        messages: list[str] = []
        for _ in range(C.CALIBRATION_MAX_ITER):
            worst: tuple[float, int, float, bool, tuple[Optional[float], Optional[float]]] = (
                0.0,
                -1,
                0.0,
                False,
                (None, None),
            )
            for node, mn, mx in self.bounds:
                i = self.pos[node.node_id]
                if mn is not None and times[i] < mn:
                    severity = (mn - times[i]) / mn
                    if severity > worst[0]:
                        worst = (severity, i, mn, True, (mn, mx))  # too young
                if mx is not None and times[i] > mx:
                    severity = (times[i] - mx) / mx
                    if severity > worst[0]:
                        worst = (severity, i, mx, False, (mn, mx))  # too old
            if worst[1] < 0 or worst[0] < 1e-12:
                return times, messages
            severity, idx, target, too_young, (row_min, row_max) = worst
            node = self.internal_nodes[idx]
            # scale the node's stem lineage so the node reaches the bound
            p = node.parent
            if p is None:
                raise ValueError(
                    f"conflicting calibration constraints: the root node cannot "
                    f"be made {'older' if too_young else 'younger'} to satisfy "
                    f"bound {target:.4g}; relax the calibrations"
                )
            if p.is_root():
                # root-level violation: adjust the global scale for this node
                current = times[idx]
                times *= target / max(current, 1e-300)
                messages.append(
                    f"scaled all times by {target / max(current, 1e-300):.4g} to "
                    f"satisfy node {node.node_id}"
                )
                continue
            parent_time = times[self.pos[p.node_id]]
            dt_stem = parent_time - times[idx]
            need = parent_time - target
            if need <= 0 or dt_stem <= 0:
                raise ValueError(
                    f"conflicting calibration constraints: node {node.node_id} "
                    f"(bound [min={_fmt_bound(row_min)}, max={_fmt_bound(row_max)}]) "
                    f"cannot be placed at {target:.6g}, because its ancestor "
                    f"node {p.node_id} sits at {parent_time:.6g} and a descendant "
                    "can never be older than its ancestor. Relax the bound on "
                    "node "
                    f"{node.node_id} or move the bound on node {p.node_id} "
                    "toward the root."
                )
            q = dt_stem / need  # rate scale (<1 slows -> older; >1 -> younger)
            self._scale_subtree(idx, q, times)
            messages.append(
                f"adjusted lineage of node {node.node_id} by rate factor {q:.4g}"
            )
        raise ValueError(
            "calibration constraints could not be satisfied within "
            f"{C.CALIBRATION_MAX_ITER} iterations; conflicting constraints at: "
            + "; ".join(
                f"node {n.node_id} [{mn}, {mx}]" for n, mn, mx in self.bounds
            )
        )

    def _check_time_order(self, times: np.ndarray) -> None:
        """Raise unless the converged ages are non-negative and nested.

        ``calibrate`` emits every tip at age 0.0, so two conditions over the
        internal vector are exactly monotonicity and non-negativity over *all*
        emitted node ages: an internal age below 0 inverts the parent--child
        pair whose child is a tip, and an internal child older than its
        internal parent inverts a nested pair.

        The violation loop needs both.  It tests the worst violator of the
        round it is fixing, so a repair that drags one lineage onto its bound
        can leave an already-inverted subtree untouched - the loop exits as soon
        as no *calibrated row* is violated any more, not as soon as the ages are
        physically ordered.  The shipped ``example_calibrations.tsv`` under
        ``method="effective"`` did exactly that: node 307 came back at
        ``-3998.56`` Ma, node 306 older than its own ancestor 305, four
        parent--child inversions and four ``NaN`` implied rates, reported as
        warnings only.  Returning such a vector would make every downstream
        number (implied rates, CIs, the written timetree) quietly wrong, so the
        solution is rejected here and the offending nodes are named.
        """
        negative = [
            self.internal_nodes[i].node_id
            for i in range(times.size)
            if times[i] < -_MONOTONIC_TOL
        ]
        inversions: list[tuple[int, int]] = []
        for i, node in enumerate(self.internal_nodes):
            parent = node.parent
            if parent is None:
                continue
            pi = self.pos.get(parent.node_id)
            if pi is None:  # pragma: no cover - ingroup parents are internal
                continue
            if times[i] > times[pi] + _MONOTONIC_TOL * max(abs(times[pi]), 1.0):
                inversions.append((node.node_id, parent.node_id))
        if not negative and not inversions:
            return
        details: list[str] = []
        if negative:
            details.append(
                f"{len(negative)} node(s) older than the present (negative age): "
                f"{sorted(negative)}"
            )
        if inversions:
            shown = "; ".join(
                f"node {child} at {times[self.pos[child]]:.6g} sits below its "
                f"ancestor node {ancestor} at {times[self.pos[ancestor]]:.6g}"
                for child, ancestor in inversions[:5]
            )
            more = f" (and {len(inversions) - 5} more)" if len(inversions) > 5 else ""
            details.append(
                f"{len(inversions)} parent-child inversion(s): {shown}{more}"
            )
        raise ValueError(
            "conflicting calibration constraints: the repair loop converged on "
            "node ages that are not nested in the tree -- " + "; ".join(details)
            + ".  This means the calibration set is unsatisfiable by a "
            "rate rescaling of this tree; relax the bounds or move them toward "
            "the root."
        )

    def _scale_subtree(self, idx: int, q: float, times: np.ndarray) -> None:
        """Scale all durations within the lineage of ``idx`` by ``1/q``."""
        node = self.internal_nodes[idx]
        p = node.parent
        assert p is not None and not p.is_root() and p.parent is not None
        parent_time = times[self.pos[p.node_id]]
        # stem: new t = parent_time - (parent_time - t)/q
        times[idx] = parent_time - (parent_time - times[idx]) / q
        stack = [idx]
        while stack:
            cur = stack.pop()
            cur_node = self.internal_nodes[cur]
            for child in cur_node.children:
                if child.is_tip():
                    continue
                ci = self.pos[child.node_id]
                times[ci] = times[cur] - (times[cur] - times[ci]) / q
                stack.append(ci)


@dataclass
class CalibrationInput:
    """Resolved calibrations attached to tree nodes (kept for API clarity)."""

    rel_result: Any
    entries: list[tuple[PhyloNode, Calibration]]


def _resolve_targets(
    tree: PhyloNode,
    calibrations: list[Calibration],
    taxon_table: Optional[dict[str, list[str]]] = None,
) -> list[tuple[PhyloNode, Calibration]]:
    id_to_node = {n.node_id: n for n in tree.walk()}
    resolved: list[tuple[PhyloNode, Calibration]] = []
    for cal in calibrations:
        if cal.node_id is not None:
            node = id_to_node.get(cal.node_id)
            if node is None or node.is_tip():
                raise KeyError(
                    f"calibration node_id {cal.node_id} is not an internal node "
                    f"of the ingroup tree"
                )
            resolved.append((node, cal))
        else:
            assert cal.taxon_set is not None
            tokens = sorted(cal.taxon_set)
            labels = set(tree.tip_labels())
            if not set(tokens) <= labels:
                # some tokens are not tips: try group-name expansion
                tokens = resolve_taxa_tokens(
                    tree, tokens, taxon_table=taxon_table, context="calibration"
                )
            result = warn_if_not_monophyletic(
                tree, tokens, name="|".join(sorted(cal.taxon_set))
            )
            node = tree.mrca(result.tips)
            if node.is_tip():
                raise KeyError(
                    f"calibration taxon set {sorted(cal.taxon_set)!r} resolves to a tip"
                )
            if set(result.tips) != set(cal.taxon_set):
                cal = Calibration(
                    taxon_set=frozenset(result.tips),
                    node_id=cal.node_id,
                    min_bound=cal.min_bound,
                    max_bound=cal.max_bound,
                    density=cal.density,
                    density_params=cal.density_params,
                )
            resolved.append((node, cal))
    return resolved


def _bounds_of(cal: Calibration) -> tuple[Optional[float], Optional[float]]:
    return cal.min_bound, cal.max_bound


def calibrate(
    times,
    calibrations: list[Calibration],
    *,
    method: Literal["bounds", "effective"] = "bounds",
    n_effective: int = C.EFFECTIVE_BOUNDS_REPLICATES,
    seed: Optional[int] = None,
    taxon_table: Optional[dict[str, list[str]]] = None,
    progress: Optional[Any] = None,
    **params: Any,
) -> CalibratedResult:
    """Convert relative times to absolute times using calibrations.

    Parameters
    ----------
    times:
        :class:`~openreltime.report.TimeResult` from :func:`openreltime.rrf_times`.
    calibrations:
        List of :class:`Calibration` (or ``parse_calibrations`` output).  A
        ``taxon_set`` may contain group names (resolved through
        ``taxon_table`` or against prefixes auto-detected from the tip
        labels); a warning is logged whenever a taxon set is not
        monophyletic.  Must not be empty, and every row must carry at least one
        finite bound (a density alone is only accepted by ``method="effective"``).
    method:
        ``"bounds"`` uses the min/max constraints directly; ``"effective"``
        derives effective bounds from densities with ``n_effective``
        replicates (msz236), also returning the empirical node-age
        distribution quantiles.  ``n_effective`` must be at least 2, since the
        effective bounds are percentiles of the replicate distribution.
    seed:
        Seed for the effective-bounds sampling.
    taxon_table:
        Optional ``{group_name: [tip, ...]}`` mapping used to expand
        group-name ``taxon_set`` tokens.
    progress:
        Optional callback ``(done, total) -> bool`` invoked after every
        effective-bounds replicate; return ``False`` to cancel, which raises
        :class:`CalibrationCancelled`.
    **params:
        Provenance recorded verbatim into ``CalibratedResult.params`` (the CLI
        and other callers pass ``tree_file``, ``calibrations_file``,
        ``outgroup``, ``input_fmt``, ``resolve_polytomy`` and
        ``outgroup_check`` through here).

    Returns
    -------
    CalibratedResult
        Absolute-time tree, the global factor ``f`` and diagnostics.  Two
        different rate vectors are reported: ``rates`` are the RRF *relative*
        lineage rates (msz236 eqs 1-4, the same quantity as ``TimeResult.rates``
        and the one :func:`openreltime.ci` uses for ``Vobs(R)``), while
        ``implied_rates`` are absolute rates per unit time obtained from the
        bare stem length divided by the calibrated duration; a lineage whose
        calibrated duration is not positive gets ``NaN`` there and is counted in
        ``params["n_lineages_dropped"]`` instead of being reported as a rate of
        0.

    Raises
    ------
    ValueError
        For an empty calibration list, a density-only row under
        ``method="bounds"``, ``n_effective < 2``, an unknown or missing
        ``density_params`` key (:data:`DENSITY_PARAM_KEYS`), and - the four
        ways a constraint set can be *unsatisfiable* - the root being asked to
        become younger/older than the tree allows, a descendant bound whose
        ancestor is already younger than that bound asks it to be, the
        violation loop exhausting :data:`openreltime._constants.CALIBRATION_MAX_ITER`
        rounds, and the converged vector failing to be nested in the tree
        (a negative age or a descendant older than its ancestor, named in the
        message).  The two nested cases raise rather than warning: a returned
        result would have violated the caller's own ``min_bound``/``max_bound``
        and silently abandoned every remaining row.
    """
    if not calibrations:
        raise ValueError(
            "no calibrations given: pass at least one Calibration (use "
            "rrf_times directly for unconstrained relative times)"
        )
    tree = times.tree
    resolved = _resolve_targets(tree, calibrations, taxon_table)

    def bounds_list(pairs):
        return [(node, *_bounds_of(cal)) for node, cal in pairs]

    warnings: list[str] = []
    # bounds actually applied to the final point estimate: the CI layer must
    # differentiate the estimator that produced the ages it qualifies, not the
    # raw rows the user typed
    used_bounds: dict[int, tuple[Optional[float], Optional[float]]] = {
        node.node_id: (cal.min_bound, cal.max_bound) for node, cal in resolved
    }
    n_rep: Optional[int] = None
    if method == "bounds":
        density_only = [
            (node, cal) for node, cal in resolved if cal.is_density_only
        ]
        if density_only:
            raise ValueError(
                f"method='bounds' cannot use density-only calibrations: node(s) "
                f"{sorted(node.node_id for node, _ in density_only)} carry a "
                "probability density but no min_bound/max_bound, and the bounds "
                "solver has no way to turn a density into a time scale (it would "
                "leave the global factor unbounded).  Either add a finite bound to "
                "those rows or re-run with method='effective', which derives "
                "msz236 effective bounds from the densities."
            )
        solver = _CalibratedTree(tree, times, bounds_list(resolved))
        f, abs_times, log = solver.solve()
        warnings.extend(log)
        effective_bounds = None
        empirical = None
    elif method == "effective":
        rng = np.random.default_rng(seed)
        dense = [(n, c) for n, c in resolved if c.density is not None]
        fixed = [(n, c) for n, c in resolved if c.density is None]
        if not dense:
            raise ValueError(
                "effective method requires at least one calibration with a "
                "density; use method='bounds' for pure min/max constraints"
            )
        records: list[dict[str, float]] = []
        keep = {n.node_id: [] for n, _ in resolved}
        # msz236 effective bounds are *percentiles of a replicate distribution*:
        # with fewer than two replicates the 2.5/97.5 pair is undefined.  Running
        # two and reporting one (the old ``max(n_effective, 2)``) made the JSON
        # record a replicate count that was never the number of solves performed.
        n_rep = int(n_effective)
        if n_rep < 2:
            raise ValueError(
                "method='effective' needs n_effective >= 2 replicates, got "
                f"{n_effective!r}: the effective bounds are the 2.5/97.5 "
                "percentiles of the replicate node-age distributions, which are "
                "undefined for a single replicate (they would collapse to the "
                "point estimate and report a zero-width interval). Use "
                "method='bounds' for a single deterministic solve."
            )
        all_nodes = [n for n in tree.walk() if not n.is_tip()]
        rep_frame = {n.node_id: [] for n in all_nodes}
        for rep in range(n_rep):
            pairs = list(fixed)
            for node, cal in dense:
                mn, mx = _sample_density_pair(cal, rng)
                # the pair is used as a hard (min, max) constraint, so both
                # dates must be positive ages (a normal density can sample
                # below zero)
                if mx <= 0:
                    raise ValueError(
                        f"the {cal.density} density at node {node.node_id} sampled "
                        f"two non-positive dates ({mn:.4g}, {mx:.4g}); calibration "
                        "ages must be positive - shift the density (offset/mean)"
                    )
                pairs.append((node, Calibration(
                    taxon_set=cal.taxon_set,
                    node_id=cal.node_id,
                    min_bound=max(mn, 0.0),
                    max_bound=mx,
                )))
            solver = _CalibratedTree(tree, times, bounds_list(pairs))
            _, abs_times, _log = solver.solve()
            for node, _cal in resolved:
                keep[node.node_id].append(float(abs_times[solver.pos[node.node_id]]))
            for n in all_nodes:
                rep_frame[n.node_id].append(float(abs_times[solver.pos[n.node_id]]))
            if progress is not None:
                if not progress(rep + 1, n_rep):
                    raise CalibrationCancelled(
                        "effective-bounds calibration cancelled by the caller "
                        f"after {rep + 1}/{n_rep} replicates"
                    )
        # effective bounds = 2.5/97.5 percentiles of calibrated-node ages
        rows = []
        for node, cal in resolved:
            values = np.array(keep[node.node_id])
            rows.append(
                {
                    "node_id": node.node_id,
                    "target": sorted(cal.taxon_set or [])
                    if cal.taxon_set
                    else f"node{node.node_id}",
                    "effective_min": float(np.quantile(values, 0.025)),
                    "effective_max": float(np.quantile(values, 0.975)),
                }
            )
        effective_bounds = pd.DataFrame(rows)
        # empirical CI over all nodes
        eframes = []
        for n in tree.walk():
            if n.is_tip():
                continue
            values = np.array(rep_frame[n.node_id])
            eframes.append(
                {
                    "node_id": n.node_id,
                    "time": float(np.median(values)),
                    "lower": float(np.quantile(values, 0.025)),
                    "upper": float(np.quantile(values, 0.975)),
                }
            )
        empirical = pd.DataFrame(eframes)
        # final point estimate with effective bounds
        eff_pairs: list[tuple[PhyloNode, Calibration]] = []
        for node, cal in resolved:
            row = effective_bounds[effective_bounds["node_id"] == node.node_id].iloc[0]
            eff_pairs.append(
                (
                    node,
                    Calibration(
                        taxon_set=cal.taxon_set,
                        node_id=cal.node_id,
                        min_bound=float(row["effective_min"]),
                        max_bound=float(row["effective_max"]),
                    ),
                )
            )
        solver = _CalibratedTree(tree, times, bounds_list(eff_pairs))
        f, abs_times, log = solver.solve()
        warnings.extend(log)
        for node, cal in eff_pairs:
            used_bounds[node.node_id] = (cal.min_bound, cal.max_bound)
    else:
        raise ValueError(f"unknown calibration method: {method!r}")

    abs_map: dict[int, float] = {tree.node_id: float(abs_times[solver.root_pos])}
    for n in tree.walk():
        if n.is_tip():
            abs_map[n.node_id] = 0.0
        else:
            abs_map[n.node_id] = float(abs_times[solver.pos[n.node_id]])

    # Two different rate vectors are reported.
    # ``rrf_rates`` are the RRF relative lineage rates, i.e. the r_j of msz236
    # eqs (1)-(4): they are defined on the folded subtree lengths and are exactly
    # the quantity rrf_times/rrf_rates report as ``Rate``.  The variance
    # decomposition of openreltime.ci (eqs 9-14) is stated for these rates.
    rrf_rates: dict[int, float] = {}
    for n in tree.walk():
        value = times.rates.get(n.node_id)
        rrf_rates[n.node_id] = float(value) if value is not None else float("nan")

    # ``implied_rates`` are absolute rates per unit calibrated time: bare stem
    # length / (parent age - node age).  A non-positive duration (zero-length
    # stem, or a node time replaced by an ancestral one by the rate-ratio guard)
    # carries no rate information at all, so it is reported as missing (NaN) and
    # counted, never as a rate of 0.0 that would silently enter the variance.
    implied_rates: dict[int, float] = {}
    n_lineages_dropped = 0
    for n in tree.walk():
        if n.is_root():
            implied_rates[n.node_id] = 1.0
            continue
        parent = n.parent
        duration = abs_map[parent.node_id] - abs_map[n.node_id]
        if duration > 0 and math.isfinite(duration):
            implied_rates[n.node_id] = float((n.blen or 0.0) / duration)
        else:
            implied_rates[n.node_id] = float("nan")
            n_lineages_dropped += 1
    if n_lineages_dropped:
        warnings.append(
            f"{n_lineages_dropped} lineage(s) have a non-positive calibrated "
            "duration; ImpliedRate is reported as missing (NaN) for them and "
            "they are excluded from rate summaries"
        )
        logger.warning(
            "%d lineage(s) have non-positive durations; their implied rate is "
            "NaN rather than 0.0 and they are dropped from rate averages",
            n_lineages_dropped,
        )

    records = [
        CalibrationRecord(
            target=sorted(c.taxon_set or [])
            if c.taxon_set
            else f"node{c.node_id}",
            node_id=n.node_id,
            min_bound=c.min_bound,
            max_bound=c.max_bound,
            density=c.density,
            density_params=dict(c.density_params),
            used_min_bound=used_bounds[n.node_id][0],
            used_max_bound=used_bounds[n.node_id][1],
        )
        for n, c in resolved
    ]
    merged = {
        "method": method,
        "n_effective": n_rep if method == "effective" else None,
        "seed": seed,
        "rate_ratio_threshold": times.params.get("rate_ratio_threshold", 20.0),
        "mean": times.params.get("mean", "geometric"),
        "n_lineages": sum(1 for n in tree.walk() if not n.is_root()),
        "n_lineages_dropped": n_lineages_dropped,
        "rate_definition": "RRF relative lineage rates (msz236 eqs 1-4)",
        "implied_rate_definition": "stem branch length / (parent age - node age)",
        **params,
    }
    return CalibratedResult(
        tree=tree,
        times=abs_map,
        rates=rrf_rates,
        time_factor=float(f),
        calibrations=records,
        effective_bounds=effective_bounds,
        empirical_ci=empirical,
        params=merged,
        warnings=warnings,
        implied_rates=implied_rates,
    )
