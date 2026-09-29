"""Plotting helpers (optional ``plot`` extra): timetrees, CIs, ddBD densities.

All functions only *save* figures via ``savefig``; no interactive display is
ever triggered.  matplotlib is an optional dependency imported lazily.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np

from .report import CalibratedResult, CIResult, DDBDResult, TimeResult
from .tree import PhyloNode

if TYPE_CHECKING:  # pragma: no cover
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

__all__ = ["plot_timetree", "plot_ci", "plot_ddbd"]

#: Diverging palette approximating RColorBrewer "RdYlBu" (10 classes)
_RATE_PALETTE = [
    "#a50026",
    "#d73027",
    "#f46d43",
    "#fdae61",
    "#fee090",
    "#e0f3f8",
    "#abd9e9",
    "#74add1",
    "#4575b4",
    "#313695",
]

#: Fallback colour for edges without a usable rate (root, zero-length lineage).
_RATE_FALLBACK = "#4575b4"


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "matplotlib is required for plotting; install it with "
            "`pip install matplotlib` (the package's [plot] extra)"
        ) from exc
    return plt


def _rate_ramp(
    rates: dict[int, float], scale: Literal["log", "linear"] = "log"
):
    """Return ``node_id -> palette index`` for the relative-rate colour ramp.

    RRF lineage rates are approximately log-distributed (a single tree spans
    0.03 to 22 in our test set), so a linear min/max normalisation piles
    almost every branch on one of the two palette extremes.  ``"log"``
    normalises on ``log10`` — the categorical equivalent of a ``LogNorm``
    colour map — which is the default; ``"linear"`` selects the plain min/max
    ramp instead.
    """
    positive = np.array([r for r in rates.values() if r > 0], dtype=float)
    if positive.size == 0:
        return {}
    if scale == "log" and positive.min() > 0 and positive.max() / positive.min() > 1.0:
        def transform(value: float) -> float:
            return math.log10(value)
    else:
        def transform(value: float) -> float:
            return value
    tvals = np.array([transform(r) for r in positive], dtype=float)
    tmin, tmax = float(tvals.min()), float(tvals.max())
    steps = len(_RATE_PALETTE)
    ramp: dict[int, int] = {}
    for node_id, value in rates.items():
        if value <= 0 or tmax <= tmin:
            ramp[node_id] = 0
            continue
        frac = (transform(value) - tmin) / (tmax - tmin)
        ramp[node_id] = min(int(np.clip(frac, 0.0, 1.0) * steps), steps - 1)
    return ramp


def plot_timetree(
    result: TimeResult | CalibratedResult,
    path: str | Path,
    *,
    use_rates: bool = True,
    rate_scale: Literal["log", "linear"] = "log",
    figsize: tuple[float, float] = (8, 6),
    dpi: int = 200,
) -> Path:
    """Draw a timetree with branches coloured by relative rates.

    Tips are ordered along the y axis; node x positions are the (relative or
    absolute) node times.  Rate colours follow :func:`_rate_ramp`
    (``rate_scale="log"`` by default).  Saved to *path* (png/pdf by
    extension).
    """
    plt = _require_matplotlib()
    tree: PhyloNode = result.tree
    times = result.times
    fig: Figure = plt.figure(figsize=figsize)
    ax: Axes = fig.add_subplot(111)

    tips = tree.tips()
    y_of: dict[int, float] = {}
    for i, tip in enumerate(tips):
        y_of[tip.node_id] = float(i)
    for node in tree.iter_postorder():
        if node.is_tip():
            continue
        kids = [y_of[c.node_id] for c in node.children]
        y_of[node.node_id] = (kids[0] + kids[-1]) / 2

    rates = getattr(result, "rates", {})
    ramp = _rate_ramp(rates, rate_scale) if use_rates else {}

    def rate_color(node: PhyloNode) -> str:
        if not use_rates or node.node_id not in ramp:
            return _RATE_FALLBACK
        return _RATE_PALETTE[ramp[node.node_id]]

    for node in tree.walk():
        if node.is_root():
            continue
        parent = node.parent
        assert parent is not None
        x0, x1 = times[parent.node_id], times[node.node_id]
        y0, y1 = y_of[parent.node_id], y_of[node.node_id]
        ax.plot([x0, x1], [y0, y0], color=rate_color(node), lw=1.0)
        if y0 != y1:
            ax.plot([x1, x1], [y0, y1], color=rate_color(node), lw=1.0)

    ax.set_yticks([y_of[t.node_id] for t in tips])
    ax.set_yticklabels([t.label or "" for t in tips], fontsize=4)
    ax.set_xlabel("Time" if isinstance(result, CalibratedResult) else "Relative time")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = Path(path)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out


def plot_ci(
    ci_result: CIResult,
    path: str | Path,
    *,
    max_nodes: int = 60,
    figsize: tuple[float, float] = (8, 6),
    dpi: int = 200,
) -> Path:
    """Plot node ages with their confidence intervals (deepest nodes first)."""
    plt = _require_matplotlib()
    frame = ci_result.table.sort_values("time", ascending=False).head(max_nodes)
    fig: Figure = plt.figure(figsize=figsize)
    ax: Axes = fig.add_subplot(111)
    y = np.arange(len(frame))[::-1]
    ax.errorbar(
        frame["time"],
        y,
        xerr=[frame["time"] - frame["lower"], frame["upper"] - frame["time"]],
        fmt="o",
        ms=3,
        lw=1,
        color="#313695",
        ecolor="#74add1",
    )
    ax.set_yticks(y)
    ax.set_yticklabels(
        [str(v) for v in frame["label"]], fontsize=5
    )
    ax.set_xlabel("Divergence time")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = Path(path)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out


def plot_ddbd(
    ddbd_result: DDBDResult,
    times: TimeResult,
    path: str | Path,
    *,
    dpi: int = 200,
    figsize: tuple[float, float] = (6, 4),
) -> Path:
    """Plot the node-time histogram with the fitted birth-death density.

    Both are drawn on the *absolute* (anchored) node-time axis, so the curve
    integrates to 1 exactly like the ``density=True`` histogram.
    """
    from .ddbd import _bd_density, _r_density

    plt = _require_matplotlib()
    internal = np.array(
        [times.times[n.node_id] for n in times.tree.walk() if not n.is_tip()]
    )
    rel = internal / internal.max()
    rel = np.clip(rel, 0, None)
    sf = ddbd_result.scale_factor
    if sf <= 0:
        raise ValueError(f"ddBD scale factor must be positive, got {sf}")
    grid, density = _r_density(rel * sf)
    fig: Figure = plt.figure(figsize=figsize)
    ax: Axes = fig.add_subplot(111)
    ax.hist(rel * sf, bins=30, density=True, color="lightgray", edgecolor="white")
    # ``_bd_density`` is a density in *relative* time u = x / sf; evaluated on
    # the plotted axis x = u * sf it must be divided by the Jacobian sf, else
    # the curve integrates to sf instead of 1 (this is what R3F's own
    # ``lines(nt*sf, bd.density.inf/sf)`` does).
    lam = ddbd_result.birth_rate * sf
    mu = ddbd_result.death_rate * sf
    fitted = (
        _bd_density(grid / sf, lam, mu, ddbd_result.sampling_frac, root_age=1.0) / sf
    )
    # non-finite density values (degenerate parameter combinations) would blow
    # up the y autoscale; matplotlib simply gaps over NaN
    ax.plot(grid, np.where(np.isfinite(fitted), fitted, np.nan), color="darkred", lw=2)
    ax.set_xlabel("Node time")
    ax.set_ylabel("Density")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = Path(path)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out
