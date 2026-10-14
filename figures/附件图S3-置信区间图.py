# -*- coding: utf-8 -*-
"""OpenRelTime Supplementary Figure S3 -- effective bounds, calibration
sensitivity and an internal inconsistency in the shipped example file.

Every number is read from the data file checked in beside this script
(``附件图S3-置信区间图-data.json``); the render path never imports
``openreltime``, so the figure reproduces without the package installed.  One
run writes a PDF exactly 178 mm wide and a 1200 dpi PNG, then checks that every
text and artwork bounding box lies inside the page rectangle.

Three panels:
  (a) one small multiple per calibration, each with its own x scale (the four
      rows span 6.56, 73.54, 33.69 and 6.65 Ma as drawn, so intervals two
      orders of magnitude apart are equally legible).  Each row states the
      prior it actually uses (family and parameters), because three of the four
      rows are not verbatim copies of ``data/examples/example_calibrations.tsv``:
      the parent density is drawn as an outline truncated to its 2.5-97.5%
      quantile interval (a quantile interval, not the density support), the
      effective bounds as a blue bar, the calibrated point as a vertical line,
      and a point clipped onto a bound is ringed.  The last row is a jointly
      feasible uniform control.
  (b) the calibration-sensitivity sweep: node age against anchored ingroup-root
      age for the eleven landmark nodes, drawn from the recorded sweep
      ``sensitivity`` (the four measured root ages 72.5, 87.5, 97.5 and
      107.5 Ma); the segments are extended to the origin because the RRF
      solution is exactly homogeneous in the root age, so the slope is the
      propagation coefficient.  The two extreme slopes are labelled.
  (c) the global time factor each row of the shipped example file implies, as a
      lollipop point plot on a log axis, with the value printed beside the point.

Run:
    python 附件图S3-置信区间图.py
    python 附件图S3-置信区间图.py --engine-root /path/to/OpenRelTime

Refresh the data (optional; needs the OpenRelTime repository and recomputes
10,000 effective-bounds replicates at the recorded seed).  The published record
is the checked-in -data.json, so the refresh writes a COPY next to it
(``附件图S3-置信区间图-data.regen.json``) unless you point it somewhere with
``--regen-out PATH`` or consciously overwrite with ``--force``:

    python 附件图S3-置信区间图.py --regen-data --engine-root /path/to/OpenRelTime
    python 附件图S3-置信区间图.py --regen-data --regen-out /tmp/s3-data.json
    python 附件图S3-置信区间图.py --regen-data --force      # replaces -data.json

--force is meant only for a deliberate re-measurement on a machine that is
allowed to become the new published record; it re-keys environment strings and
re-measures everything, so it is never the way to check the printed numbers.
"""
import argparse
import json
import math
import os
import sys

import numpy as np
from matplotlib.transforms import Bbox

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
PROJECT_ROOT = os.path.dirname(HERE)

STEM = "附件图S3-置信区间图"
DATA_JSON = STEM + "-data.json"
# --regen-data writes here unless --force is given: the checked-in DATA_JSON is
# the published record and must not be clobbered by a refresh on some other
# machine (review item S2-i/F3-7, applied here to S3 as well).
REGEN_JSON = STEM + "-data.regen.json"
ENV_ENGINE = "OPENRELTIME_ENGINE_ROOT"
# what --regen-data needs out of the repository, and the recorded settings the
# shipped numbers were produced with
EXAMPLE_TREE = os.path.join("data", "examples", "example.nwk")
EXAMPLE_CALIBS = os.path.join("data", "examples", "example_calibrations.tsv")
OUTGROUP = ["Ornithorhynchus_anatinus", "Zaglossus_bruijni", "Tachyglossus_aculeatus"]
SEED = 20260906
N_EFFECTIVE = 10000
NORMALIZE = True
CROWN_SCENARIOS = [{"mid_ma": 72.5, "tag": "65-80"}, {"mid_ma": 87.5, "tag": "80-95"},
                   {"mid_ma": 97.5, "tag": "90-105"}, {"mid_ma": 107.5, "tag": "100-115"}]
ANCHOR_MA = 97.5                                  # the scenario panel (b) is read at
# panel (a): (row label, tip pair addressed by taxon_set, density, density_params,
# extra keys).  The two rows marked bound-only carry a density here that the
# shipped TSV does not have; the fourth row is the jointly feasible control.
ROWS = [
    ("Homo-Pan", ("Homo_sapiens", "Pan_troglodytes"), "exponential",
     {"offset": 5.0, "mean": 1.5}, "bound-only row of the TSV, given a density here"),
    ("Elephas-Loxodonta", ("Elephas_maximus", "Loxodonta_africana"), "exponential",
     {"offset": 60.0, "mean": 20.0}, "density row of the shipped TSV"),
    ("Dasypus-Choloepus", ("Dasypus_novemcinctus", "Choloepus_didactylus"), "uniform",
     {"min": 70.0, "max": 105.0}, "bound-only row of the TSV, read as a uniform density"),
    ("Homo-Pan (uniform control)", ("Homo_sapiens", "Pan_troglodytes"), "uniform",
     {"min": 5.0, "max": 12.0}, "control, jointly feasible with itself"),
]
# the eleven landmark nodes of panel (b), addressed by the taxon set whose MRCA
# they are; the label is the caption's own "Genus species - Genus species" name
# and is NOT derived from the node's tip set (only 2 of the 11 MRCAs are
# two-tip nodes -- Homo-Pan has 3 tips, Sus-Bos 63 -- so the MRCA of the named
# pair is the right way to address them)
LANDMARKS = [
    ("Homo sapiens - Pan troglodytes", ("Homo_sapiens", "Pan_troglodytes")),
    ("Homo sapiens - Gorilla gorilla", ("Homo_sapiens", "Gorilla_gorilla")),
    ("Homo sapiens - Pongo abelii", ("Homo_sapiens", "Pongo_abelii")),
    ("Homo sapiens - Mus musculus", ("Homo_sapiens", "Mus_musculus")),
    ("Mus musculus - Rattus norvegicus", ("Mus_musculus", "Rattus_norvegicus")),
    ("Loxodonta africana - Elephas maximus", ("Loxodonta_africana", "Elephas_maximus")),
    ("Dasypus novemcinctus - Choloepus didactylus",
     ("Dasypus_novemcinctus", "Choloepus_didactylus")),
    ("Ursus maritimus - Ursus arctos", ("Ursus_maritimus", "Ursus_arctos")),
    ("Balaenoptera musculus - Megaptera novaeangliae",
     ("Balaenoptera_musculus", "Megaptera_novaeangliae")),
    ("Tursiops truncatus - Delphinus capensis",
     ("Tursiops_truncatus", "Delphinus_capensis")),
    ("Sus scrofa - Bos taurus", ("Sus_scrofa", "Bos_taurus")),
]


def locate_engine_root(explicit=None):
    """Find the OpenRelTime repository by content, never by folder name.

    The repository supplies ``data/examples/example.nwk`` and the shipped
    calibration table to ``--regen-data``.  As shipped to a journal the sibling
    folders may be renamed or absent, so every candidate is probed for the
    package directory plus the example tree rather than matched by name.
    """
    for given in (explicit, os.environ.get(ENV_ENGINE), PROJECT_ROOT, os.getcwd()):
        if not given:
            continue
        path = os.path.abspath(os.path.expanduser(given))
        if os.path.isdir(os.path.join(path, "openreltime")) and \
                os.path.isfile(os.path.join(path, EXAMPLE_TREE)):
            return path
    try:
        siblings = sorted(os.listdir(PROJECT_ROOT))
    except OSError:
        siblings = []
    for name in siblings:
        path = os.path.join(PROJECT_ROOT, name)
        if os.path.isdir(os.path.join(path, "openreltime")) and \
                os.path.isfile(os.path.join(path, EXAMPLE_TREE)):
            return path
    return None

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter, LogLocator, MaxNLocator, NullFormatter  # noqa: E402
from scipy.stats import expon as _expon, uniform as _uniform  # noqa: E402

# --- shared style constants, identical to the S2 script -------------------- #
MM = 1.0 / 25.4
COLUMN_MM = 178.0
FS_PANEL = 9.0
FS_AXIS = 8.5
FS_TICK = 8.0
FS_SERIES = 7.0
PANEL_XY = (-0.02, 1.02)
# Okabe-Ito / colour-blind-safe palette
BLUE = "#0072B2"
VERM = "#D55E00"
GREEN = "#009E73"
ORANGE = "#E69F00"
GREY = "#555555"
DENSITY_FILL = "#dcdcdc"
DENSITY_EDGE = "#9a9a9a"
MARK = "#1b3a56"

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
    "font.size": FS_TICK, "axes.labelsize": FS_AXIS,
    "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
    "axes.linewidth": 0.8, "axes.edgecolor": "#333333",
    "savefig.dpi": 1200, "figure.dpi": 100,
})

# one canonical short name per landmark node (review F4-9); en dashes, as the
# caption spells them, and "Lox-Ele" replaced by the caption's order
ABBR = {"Homo sapiens - Pan troglodytes": "Homo\u2013Pan",
        "Homo sapiens - Gorilla gorilla": "Homo\u2013Gorilla",
        "Homo sapiens - Pongo abelii": "Homo\u2013Pongo",
        "Homo sapiens - Mus musculus": "Homo\u2013Mus",
        "Mus musculus - Rattus norvegicus": "Mus\u2013Rattus",
        "Loxodonta africana - Elephas maximus": "Lox\u2013Ele",
        "Dasypus novemcinctus - Choloepus didactylus": "Das\u2013Cho",
        "Ursus maritimus - Ursus arctos": "Urs\u2013Urs",
        "Balaenoptera musculus - Megaptera novaeangliae": "Bal\u2013Meg",
        "Tursiops truncatus - Delphinus capensis": "Tur\u2013Del",
        "Sus scrofa - Bos taurus": "Sus\u2013Bos"}
# family and parameters of the prior behind each row of panel (a), as reported
# by the data file (density / density_params of effective_bounds)
PRIOR_NOTE = {
    "Homo-Pan": "prior: exponential\noffset 5 Ma, mean 1.5 Ma",
    "Elephas-Loxodonta": "prior: exponential\noffset 60 Ma, mean 20 Ma",
    "Dasypus-Choloepus": "prior: uniform\n70\u2013105 Ma",
    "Homo-Pan (uniform control)": "prior: uniform\n5\u201312 Ma (control)",
}
# review F6-5: 6.0 pt was below the engine's own 7.5 pt floor (make_figures.py)
FS_PRIOR = 7.0

ROW_LABEL = {"Homo-Pan": "Homo\u2013Pan (exp)",
             "Elephas-Loxodonta": "Lox\u2013Ele (exp)",
             "Dasypus-Choloepus": "Dasypus\u2013Chol. (U)",
             "Homo-Pan (uniform control)": "Homo\u2013Pan ctrl (U)"}
# panel (c) short names, keyed by the data-file row name, same order as ROW_LABEL
C_ROW_LABEL = {"Homo-Pan": "Homo\u2013Pan", "Elephas-Loxodonta": "Lox\u2013Ele",
               "Dasypus-Choloepus": "Dasypus\u2013Chol."}
# F4-8: panel (b) plots the recorded sweep, so the in-panel note has to say which
# root ages were actually run and why the segments still reach (0, 0).  Wrapped
# to <= 24 characters so the block stays clear of the steepest line.
SWEEP_NOTE = ("measured root ages\n72.5, 87.5, 97.5, 107.5\n"
              "Ma; lines extended to\n0 Ma: RRF is exactly\n"
              "homogeneous in root age")
SWEEP_HEADROOM = 1.40          # y headroom over the deepest sweep value, so the
                               # note above has a clear band (was 1.10)



# caption label of panel (b) and data-file row name of panel (c) for a taxon
# pair.  --regen-data keys the data file by these strings: a landmark is
# addressed by the MRCA of the named pair, never by the pair being its node's
# complete tip set.  ROWS[0] and ROWS[3] share the Homo-Pan pair, so the control
# row is excluded from ROW_BY_PAIR (panel (c) reports the three shipped rows).
LABEL_BY_PAIR = {frozenset(pair): label for label, pair in LANDMARKS}
ROW_BY_PAIR = {frozenset(pair): name for name, pair, _d, _p, _w in ROWS
               if "control" not in name}


# --------------------------------------------------------------------------- #
# Only --regen-data imports openreltime and recomputes the data file.  The
# default render path never enters this function, which is what keeps the figure
# reproducible without the package installed (mirrors the S2 script).
#
# Review P0-5 (BLOCKER): this branch could not run at all, for two reasons:
#   (1) `import openreltime` stood ABOVE the locate_engine_root()/sys.path
#       block, so a source checkout given by --engine-root was never on the
#       import path: ModuleNotFoundError.  Fixed below, exactly as the S2
#       script had already fixed it.
#   (2) the landmark keys were derived from each node's own tip set and the
#       derivation raised "RuntimeError: expected a two-tip landmark node,
#       got 3 tips".  Only 2 of the 11 landmark MRCAs are two-tip nodes
#       (Homo-Pan takes in Pan_paniscus, Sus-Bos takes in 63 tips), so the very
#       first landmark aborted the run and nothing was ever written.  A
#       landmark is a label -- the MRCA of the named pair -- so the label from
#       LANDMARKS/ROWS is now used directly; the derived key also mis-ordered
#       the two genuinely two-tip nodes ("Elephas - Loxodonta" by tree order
#       against the caption's "Loxodonta - Elephas").
# The refresh no longer replaces the published record: it writes REGEN_JSON, or
# the --regen-out PATH, and only --force writes back into -data.json.
# --------------------------------------------------------------------------- #
def regen_data(engine_root=None, out_json=None, force=False):  # pragma: no cover
    """Recompute every number of the data file by calling the engine.

    Panel (a) comes from one joint run of the three density rows of the variant
    calibration set, ``calibrate(..., method='effective', n_effective=10000,
    seed=20260906)``, plus a separate run of the jointly feasible control row;
    the sweep of panel (b) is the normalised relative RRF times (root = 1)
    rescaled by each of the four root ages, cross-checked against four genuine
    single-point root calibrations; the factors of panel (c) come from the rows
    of the shipped example table.
    """
    # Resolve the repository FIRST and put it on sys.path: --engine-root points
    # at a source checkout, not an installed distribution, so importing
    # openreltime before the lookup failed even when the flag was correct.
    root = locate_engine_root(engine_root)
    if root is None:
        raise SystemExit(
            "error: --regen-data needs the OpenRelTime repository (it reads\n"
            "  %s and %s).  Pass --engine-root PATH, or set %s."
            % (EXAMPLE_TREE, EXAMPLE_CALIBS, ENV_ENGINE))
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    import openreltime as ort
    from openreltime.calibrate import Calibration, density_quantiles
    out_json = out_json or (DATA_JSON if force else REGEN_JSON)
    tree = ort.read_tree(os.path.join(root, EXAMPLE_TREE), outgroup=OUTGROUP)
    # NORMALIZE is the setting the record states (note_on_disclosure): the root
    # time is scaled to 1, which is what relative_times_normalized holds.
    times = ort.rrf_times(tree, normalize=NORMALIZE)
    root_id = tree.node_id                       # the ingroup (Theria) crown
    # with normalize=True ``times.times`` are the relative RRF times scaled so
    # that the root is 1; those ARE panel (b)'s relative times, and the anchored
    # runs below only confirm that a single root calibration rescales them

    # ---- panel (a): the effective-bounds machinery ------------------------- #
    # The checked-in record is a JOINT run of the three density rows plus a
    # separate run of the jointly feasible control row.  Four independent
    # single-row runs -- what this function used to do -- cannot reproduce it:
    # alone, the Homo-Pan exponential gives effective_max 9.26 Ma rather than
    # the recorded 11.599 Ma, and no point lands on a bound, so the three
    # clipped-to-bound flags of panel (a) were uncheckable.  Jointly, the
    # conflicting bounds drive the repair loop, which is the point of the panel.
    effective_bounds, support = [], {}
    groups = [[r for r in ROWS if "control" not in r[0]],
              [r for r in ROWS if "control" in r[0]]]
    measured = {}
    for group in groups:
        if not group:
            continue
        res = ort.calibrate(
            times, [Calibration(taxon_set=frozenset(p), density=d,
                                density_params=dict(pr))
                    for _n, p, d, pr, _w in group],
            method="effective", n_effective=N_EFFECTIVE, seed=SEED)
        by_node = {int(r["node_id"]): r for _i, r in res.effective_bounds.iterrows()}
        for name, pair, _d, _p, _w in group:
            measured[name] = (res, by_node[tree.mrca(list(pair)).node_id])
    for name, pair, dens, params, _why in ROWS:
        res, row = measured[name]
        node_id = int(row["node_id"])
        lo, hi = float(row["effective_min"]), float(row["effective_max"])
        point = float(res.times[node_id])
        tol = 1e-9 * max(1.0, abs(hi))
        on_lo, on_hi = abs(point - lo) <= tol, abs(point - hi) <= tol
        inside = lo - tol <= point <= hi + tol
        if "control" in name:
            # the control row of the checked-in file states its density before
            # the bounds, so the refresh keeps that key order and diffs clean
            qlo, qhi = density_quantiles(dens, params)
            rec = {"node": name, "node_id": node_id, "density": dens,
                   "density_params": dict(params),
                   "density_quantile_2p5_97p5": [qlo, qhi],
                   "eff_min": lo, "eff_max": hi, "point": point,
                   "on_lower_bound": on_lo, "on_upper_bound": on_hi,
                   "inside_effective": inside, "strictly_interior": lo < point < hi,
                   "role": "control_jointly_feasible"}
        else:
            rec = {"node": name, "node_id": node_id, "eff_min": lo, "eff_max": hi,
                   "point": point, "on_lower_bound": on_lo, "on_upper_bound": on_hi,
                   "inside_effective": inside,
                   "density": dens, "density_params": dict(params)}
        effective_bounds.append(rec)
        support[name] = list(density_quantiles(dens, params))

    # ---- panel (b): the sensitivity sweep, anchored at the ingroup root ---- #
    landmark_ids = {}
    for label, pair in LANDMARKS:
        # address the landmark by the caption's label, i.e. by the MRCA of the
        # named pair; the MRCA may subtend any number of tips
        node = tree.mrca(list(pair))
        landmark_ids[label] = node.node_id
        print("  landmark %-46s MRCA over %d tips" % (label, len(node.tips())))
    # the normalised relative RRF times (root = 1, normalize as recorded), which
    # is what the data file's relative_times_normalized key holds
    relative_times = {label: float(times.times[nid]) for label, nid in landmark_ids.items()}
    relative_times["Theria crown (ingroup root)"] = float(times.times[root_id])
    # the sweep of a single root calibration is that normalisation rescaled by
    # the anchor, exactly: age = relative_time x root_age (Supplementary Text S4)
    sensitivity = {}
    for scen in CROWN_SCENARIOS:
        anchor = scen["mid_ma"]
        sensitivity[scen["tag"]] = dict(
            {label: relative_times[label] * anchor for label in landmark_ids},
            **{"Theria crown (ingroup root)": anchor})
    # cross-check the recorded construction against four genuine anchored runs
    worst = 0.0
    for scen in CROWN_SCENARIOS:
        res = ort.calibrate(times, [Calibration(node_id=root_id, min_bound=scen["mid_ma"],
                                                max_bound=scen["mid_ma"])], method="bounds")
        for label, nid in landmark_ids.items():
            got, exp = float(res.times[nid]), sensitivity[scen["tag"]][label]
            worst = max(worst, abs(got - exp) / exp)
    print("  anchored runs agree with the recorded sweep to %.2e relative "
          "(max over 44 nodes x 4 anchors)" % worst)
    # the propagation coefficient IS the normalised relative time (the caption
    # says so); the checked-in file stores the same quantity twice and its
    # propagation_coefficients copy sits 1-4 ulp off its own
    # relative_times_normalized for 10 of the 11 nodes.  No arithmetic route on
    # the recorded engine reproduces that offset, so the refresh writes the
    # consistent value and the difference is reported, not hidden.
    propagation = {k: v for k, v in relative_times.items()
                   if k != "Theria crown (ingroup root)"}

    # ---- panel (c): the factor each shipped example row implies ------------ #
    shipped = ort.parse_calibrations(os.path.join(root, EXAMPLE_CALIBS))
    rep_age, implied, tree_implied = {}, {}, {}
    for cal in shipped:
        pair = frozenset(cal.taxon_set or ())
        if pair not in ROW_BY_PAIR:
            raise SystemExit(
                "error: the example row %s cannot be keyed by a caption label; "
                "add it to ROWS/LANDMARKS or drop --regen-data." % (sorted(pair),))
        name = ROW_BY_PAIR[pair]          # the row name the data file uses
        label = LABEL_BY_PAIR[pair]       # the panel (b) landmark key
        node = tree.mrca(sorted(pair))    # validated, not re-labelled
        if cal.density == "exponential":
            age = cal.density_params["offset"] + cal.density_params["mean"]
            rule = "exponential(offset %g, mean %g): offset+mean" % (
                cal.density_params["offset"], cal.density_params["mean"])
        else:
            age = (cal.min_bound + cal.max_bound) / 2.0
            rule = "uniform(TSV %g-%g): midpoint" % (cal.min_bound, cal.max_bound)
        rel_time = relative_times[label]
        print("  example row %-20s MRCA over %d tips" % (name, len(node.tips())))
        rep_age[name] = {"value_Ma": float(age), "rule": rule}
        implied[name] = float(age) / rel_time
        tree_implied[name] = float(rel_time * ANCHOR_MA)
    spread = max(implied.values()) / min(implied.values())

    # the prose keys (note_on_disclosure, conclusion, ...) come from the
    # checked-in record, never from this machine
    with open(DATA_JSON, "r", encoding="utf-8") as fh:
        payload = json.load(fh)
    payload.update({
        "generated_by": os.path.basename(__file__),
        "openreltime_version": getattr(ort, "__version__", "unknown"),
        "seed": SEED, "n_effective": N_EFFECTIVE, "normalize": NORMALIZE,
        "crown_scenarios": CROWN_SCENARIOS,
        "relative_times_normalized": relative_times,
        "effective_bounds": effective_bounds,
        "density_support_2p5_97p5": support,
        "sensitivity": sensitivity,
        "propagation_coefficients": propagation,
        "example_calibration_consistency": {
            "implied_global_time_factor": implied,
            "spread_ratio": float(spread),
            "tree_implied_ages_at_crown_97_5Ma": tree_implied,
            "conclusion": payload.get("example_calibration_consistency", {}).get(
                "conclusion", ""),
            "representative_age": rep_age,
        },
    })
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1, ensure_ascii=False)
    if os.path.abspath(out_json) == os.path.abspath(DATA_JSON):
        print("refreshed (forced, published record replaced)", out_json)
    else:
        print("refreshed a COPY at", out_json, "-- the published record",
              DATA_JSON, "is untouched; pass --force to replace it")


def plain_log_ticks(axis, subs=(1.0,)):
    """Label a log axis with plain decimal numbers instead of 10^n superscripts."""
    vmin, vmax = axis.get_view_interval()
    if subs == (1.0,):
        lo = int(math.floor(math.log10(max(vmin, 1e-12))))
        hi = int(math.ceil(math.log10(vmax)))
        ticks = [10.0 ** k for k in range(lo, hi + 1) if vmin <= 10.0 ** k <= vmax]
    else:
        ticks = [s * 10.0 ** k for k in range(int(math.floor(math.log10(vmin))),
                                              int(math.ceil(math.log10(vmax))) + 1)
                 for s in subs if vmin <= s * 10.0 ** k <= vmax]
    axis.set_ticks(sorted(ticks))
    axis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    axis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(0.1, 10.0001, 0.1), numticks=100))
    axis.set_minor_formatter(NullFormatter())


def density_outline(row, quantile_lo, quantile_hi, n=200):
    """(x, normalised density) of the parent prior, cut at its 2.5-97.5% interval."""
    dens = row.get("density", "uniform")
    p = row.get("density_params", {})
    x = np.linspace(quantile_lo, quantile_hi, n)
    if dens == "exponential":
        pdf = _expon.pdf(x, loc=p.get("offset", 0.0), scale=p.get("mean", 1.0))
    else:
        pdf = _uniform.pdf(x, loc=p.get("min", 0.0), scale=p.get("max", 1.0) - p.get("min", 0.0))
    mx = pdf.max()
    return x, (pdf / mx if mx > 0 else pdf)


def check_text_overlaps(fig, renderer):
    """Fail the render when any two texts, or a text and a data artist, collide.

    The figures previously passed the page-bounds check while still printing a
    statistics box on top of the points it described; page containment alone is
    not sufficient.
    """
    def boxes():
        out = []
        for t in fig.texts:
            out.append(("figtext:" + repr(t.get_text()[:20]), t.get_window_extent()))
        for i, ax in enumerate(fig.axes):
            items = list(ax.texts) + [ax.xaxis.label, ax.yaxis.label, ax.title]
            items += list(ax.get_xticklabels()) + list(ax.get_yticklabels())
            for t in items:
                if not t.get_text().strip():
                    continue
                out.append(("p%d:%s" % (i, repr(t.get_text()[:20].replace("\n", "/"))),
                            t.get_window_extent(renderer)))
            if ax.get_legend() is not None:
                out.append(("p%d:legend" % i, ax.get_legend().get_window_extent(renderer)))
        return out

    def artists():
        """Yield (panel, name, [display boxes]) for every plotted datum."""
        out = []
        for i, ax in enumerate(fig.axes):
            for col in ax.collections:
                try:                                     # only marker collections
                    offs = col.get_offsets()
                    sizes = col.get_sizes()
                except Exception:
                    continue
                if offs is None or len(offs) == 0 or len(sizes) == 0:
                    continue
                pts = ax.transData.transform(offs)
                size = 0.5 * float(np.max(sizes)) ** 0.5
                boxes = [Bbox([[x - size, y - size], [x + size, y + size]])
                         for x, y in pts]
                out.append(("p%d:scatter" % i, boxes))
            for ln in ax.lines:
                # a polyline's bounding box covers the whole panel, which would
                # flag any nearby text; test its vertices instead
                if ln.get_marker() in ("", "None", None) or len(ln.get_xdata()) == 0:
                    continue
                try:
                    xy = np.column_stack([
                        np.asarray(ln.get_xdata(), dtype=float),
                        np.asarray(ln.get_ydata(), dtype=float)])
                    pts = ax.transData.transform(xy)
                except (TypeError, ValueError):
                    continue
                ms = ln.get_markersize() * fig.dpi / 72.0
                out.append(("p%d:markers" % i,
                            [Bbox([[x - ms / 2, y - ms / 2], [x + ms / 2, y + ms / 2]])
                             for x, y in pts if np.isfinite(x) and np.isfinite(y)]))
        return out

    bad = []
    tb = boxes()
    for i in range(len(tb)):
        for j in range(i + 1, len(tb)):
            a, b = tb[i][1], tb[j][1]
            ov = min(a.x1, b.x1) - max(a.x0, b.x0)
            oy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ov > 2.0 and oy > 2.0:              # tolerate hairline touching
                bad.append(("text/text", tb[i][0], tb[j][0], round(ov, 1), round(oy, 1)))
    for name, bb in tb:
        if not name.startswith("p0:") and not name.startswith("p1:"):
            continue
        for aname, dots in artists():
            if aname.split(":")[0] != name.split(":")[0]:
                continue
            if dots is None:
                continue
            for d in dots:
                db = d.bounds if hasattr(d, "bounds") else None
                if db is None:
                    continue
                ov = min(bb.x1, db[0] + db[2]) - max(bb.x0, db[0])
                oy = min(bb.y1, db[1] + db[3]) - max(bb.y0, db[1])
                if ov > 1.5 and oy > 1.5:
                    bad.append(("text/data", name, aname, round(ov, 1), round(oy, 1)))
                    break
    return bad



def check_column_territories(fig, renderer, tol_px=1.0):
    """Fail when a panel's own text reaches into a neighbouring column.

    Overlap checking alone still lets a y-axis label or an annotation cross the
    gutter and sit above the adjacent panel's data area, which reads as one
    panel labelling another.  Each column gets a territory bounded by the
    mid-point of the gutters either side; every text owned by that column must
    fall inside it.
    """
    w_px, h_px = fig.canvas.get_width_height()
    cols = {}
    for ax in fig.axes:
        pos = ax.get_position()
        key = round(pos.x0, 3)
        for k in list(cols):
            if abs(k - key) < 0.02:
                key = k
                break
        cols.setdefault(key, []).append(ax)
    order = sorted(cols)
    edges = []
    for i, key in enumerate(order):
        left = 0.0 if i == 0 else (cols[order[i - 1]][0].get_position().x1 + key) / 2.0
        right = (1.0 if i == len(order) - 1
                 else (key + cols[order[i]][0].get_position().x1
                       + cols[order[i + 1]][0].get_position().x0) / 2.0
                       if i + 1 < len(order) else 1.0)
        if i == len(order) - 1:
            right = 1.0
        edges.append((key, left, right))
    bad = []
    for key, left, right in edges:
        x0, x1 = left * w_px, right * w_px
        for ax in cols[key]:
            items = list(ax.texts) + [ax.xaxis.label, ax.yaxis.label, ax.title]
            items += list(ax.get_xticklabels()) + list(ax.get_yticklabels())
            if ax.get_legend() is not None:
                items.append(ax.get_legend())
            for t in items:
                txt = t.get_text() if not hasattr(t, "get_lines") else "LEGEND"
                if not (txt or hasattr(t, "get_lines")):
                    continue
                if not txt and not hasattr(t, "get_lines"):
                    continue
                bb = t.get_window_extent(renderer)
                if bb.x0 < x0 - tol_px or bb.x1 > x1 + tol_px:
                    bad.append((round(key, 3), repr((txt or "")[:26].replace("\n", "/")),
                                round(bb.x0, 1), round(bb.x1, 1),
                                (round(x0, 1), round(x1, 1))))
    return bad


def main(out_dir=None):
    out_dir = HERE if out_dir is None else out_dir
    with open(DATA_JSON, "r", encoding="utf-8") as fh:  # lives beside the script
        payload = json.load(fh)
    sens = payload["sensitivity"]
    rel = payload["relative_times_normalized"]
    crown_mid = {t["tag"]: t["mid_ma"] for t in payload["crown_scenarios"]}
    rows = payload["effective_bounds"]
    dens_support = payload["density_support_2p5_97p5"]
    ecc = payload["example_calibration_consistency"]
    implied = ecc["implied_global_time_factor"]
    rep = ecc["representative_age"]
    spread = ecc["spread_ratio"]

    fig = plt.figure(figsize=(COLUMN_MM * MM, COLUMN_MM * MM * 0.60))
    # three visually equal columns (review note: the panels must read as the same
    # width); the gutter carries panel (c)'s row labels, so wspace stays generous
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.0],
                          left=0.168, right=0.988, top=0.882, bottom=0.168, wspace=0.46)
    gsa = gs[0].subgridspec(len(rows), 1, hspace=0.42)
    a_axes = [fig.add_subplot(gsa[i]) for i in range(len(rows))]
    ax_b = fig.add_subplot(gs[1])
    ax_c = fig.add_subplot(gs[2])

    # ---------------- (a) small multiples: prior outline, bounds, point --------- #
    spans = []
    for i, (ax, row) in enumerate(zip(a_axes, rows)):
        nm = row["node"]
        qlo, qhi = dens_support[nm]
        dx, dh = density_outline(row, qlo, qhi)
        # C12: the four rows are not all rows of the shipped TSV, so the prior
        # each one actually uses is stated inside its own panel.
        ax.text(0.985, 0.98, PRIOR_NOTE[row.get("prior_note", nm)], transform=ax.transAxes,
                fontsize=FS_PRIOR, ha="right", va="top", color="#333333", zorder=7,
                linespacing=1.25,
                bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="none", alpha=0.88))
        ax.fill_between(dx, -dh, dh, color=DENSITY_FILL, edgecolor=DENSITY_EDGE,
                        lw=0.6, zorder=1)
        ax.plot([row["eff_min"], row["eff_max"]], [0, 0], "-", color=BLUE, lw=4.6,
                solid_capstyle="butt", zorder=3)
        ax.plot([row["point"], row["point"]], [-0.55, 0.55], "-", color=MARK, lw=1.4, zorder=5)
        if row.get("on_lower_bound") or row.get("on_upper_bound"):
            ax.plot([row["point"]], [0], "o", ms=4.4, mfc="none", mec=MARK, mew=1.2, zorder=6)
            tag = "point clipped to bound"
        else:
            tag = "point interior (jointly feasible)"
        span_lo = min(qlo, row["eff_min"], row["point"])
        span_hi = max(qhi, row["eff_max"], row["point"])
        # F4-4: the caption justifies the independent x scales with the measured
        # spans, so they are recorded here rather than rounded to "6 Ma / 73 Ma"
        spans.append((ROW_LABEL.get(nm, nm), span_hi - span_lo))
        ax.set_xlim(span_lo * 0.97, span_hi * 1.03)   # show where the blue
        # bar pokes out of the grey outline: expected, they are different intervals
        ax.set_ylim(-0.92, 0.92)
        ax.set_yticks([0]); ax.set_yticklabels([ROW_LABEL.get(nm, nm)], fontsize=FS_SERIES)
        ax.tick_params(axis="y", length=0)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3, prune="both"))
        plain_linear(ax)
        ax.grid(axis="x", lw=0.4, color="#ededed", zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        if i == len(rows) - 1:
            ax.set_xlabel("Divergence age (Ma)", fontsize=FS_AXIS)
        else:
            ax.set_xticklabels([])

    # ---------------- (b) the calibration-sensitivity sweep --------------------- #
    # drawn from the RECORDED sweep payload["sensitivity"] -- the four measured
    # root ages of crown_scenarios x the eleven landmark nodes -- not from a
    # rescaling of the relative times (review F4-8; the two are numerically
    # identical here, but the record must be what is plotted).  The segments to
    # the origin are an extrapolation, stated in-panel: the RRF solution is
    # exactly homogeneous in the root age, so the slope is the propagation
    # coefficient and equals the normalised relative time.
    nodes = [k for k in rel if k != "Theria crown (ingroup root)"]
    tags = [t["tag"] for t in payload["crown_scenarios"]]
    xsc = [crown_mid[t] for t in tags]
    xmax = max(xsc)
    sweep = {k: [sens[t][k] for t in tags] for k in nodes}
    for k in nodes:                              # 11 lines, all through the origin
        ys = sweep[k]
        ax_b.plot([0.0] + xsc, [0.0] + ys, "-", color="#8a8a8a", lw=1.0, zorder=2)
        ax_b.plot(xsc, ys, "o", ms=3.0, color="#8a8a8a", zorder=3)
    # highlight the steepest and shallowest slopes as end-point anchors
    top = max(max(ys) for ys in sweep.values())
    hi = max(nodes, key=lambda k: rel[k]); lo = min(nodes, key=lambda k: rel[k])
    for k, col in ((hi, BLUE), (lo, ORANGE)):
        ys = sweep[k]
        ax_b.plot([0.0] + xsc, [0.0] + ys, "-", color=col, lw=1.5, zorder=4)
        ax_b.plot(xsc, ys, "o", ms=3.6, color=col, zorder=5)
        lab_y = max(ys[-1], top * 0.055)         # lift clear of the axis frame
        # slopes to five decimals: at three none of the eleven printed ages is
        # recomputable as slope x 97.5 Ma (review P1-21 / F4-5)
        ax_b.annotate("%s  %.5f" % (ABBR.get(k, k), rel[k]), xy=(xmax + 1.0, ys[-1]),
                      xytext=(xmax + 5.5, lab_y), ha="left", va="center",
                      fontsize=FS_SERIES, color=col,
                      arrowprops=(dict(arrowstyle="-", color=col, lw=0.6, shrinkA=1,
                                       shrinkB=2)
                                  if lab_y > ys[-1] * 1.5 else None))
    ax_b.plot([0], [0], marker="*", ms=9, color="black", zorder=6)
    ax_b.set_xlim(0, 152)
    # the sweep still hugs the data; the extra band above carries SWEEP_NOTE
    ax_b.set_ylim(0, top * SWEEP_HEADROOM)
    ax_b.set_xlabel("Anchored ingroup-root age (Ma)")
    ax_b.set_ylabel("")
    # F4-8: which root ages were actually run, and why the lines reach (0, 0).
    # Kept inside the empty wedge above the steepest line so that it cannot
    # touch the data or the slope labels.
    ax_b.text(0.02, 0.97, SWEEP_NOTE, transform=ax_b.transAxes, fontsize=FS_SERIES,
              ha="left", va="top", color="#333333", zorder=7, linespacing=1.3,
              bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="none", alpha=0.88))

    ax_b.grid(lw=0.4, color="#ededed", zorder=0); ax_b.set_axisbelow(True)

    # ---------------- (c) lollipop points (in place of log bars) ---------------- #
    names = ["Homo-Pan", "Elephas-Loxodonta", "Dasypus-Choloepus"]
    vals = [implied[nm] for nm in names]
    vmin, vmax = min(vals), max(vals)
    xlo, xhi = vmin / 1.15, vmax * 1.15          # data-driven padding, not a fixed range
    ax_c.set_xscale("log")
    for i, v in enumerate(vals):
        ax_c.plot([xlo, v], [i, i], "-", color=BLUE, lw=1.0, zorder=2)
        ax_c.plot([v], [i], "o", ms=6.5, color=BLUE, mec="#004d78", mew=0.6, zorder=3)
        ax_c.text(v * 1.10, i, "%.0f" % v, va="center", ha="left",
                  fontsize=FS_SERIES, color="#222222")
        ax_c.text(xlo * 1.05, i + 0.26, C_ROW_LABEL.get(names[i], names[i]),
                  va="bottom", ha="left", fontsize=FS_SERIES, color="#333333")
    ax_c.set_yticks(range(len(names)))
    ax_c.set_yticklabels([])          # the row names are printed inside the panel
    ax_c.tick_params(axis="y", length=0)
    ax_c.tick_params(axis="x", which="minor", length=2.2, width=0.5)
    ax_c.set_xlim(xlo, xhi * 2.35); ax_c.set_ylim(-0.55, len(names) - 0.25)
    ax_c.set_xlabel("")
    plain_log_ticks(ax_c.xaxis, subs=(1, 5))
    ax_c.set_xticks([200, 1000, 4000])
    ax_c.set_xlim(xlo, xhi * 1.85)
    ax_c.grid(axis="x", lw=0.4, color="#ededed", zorder=0); ax_c.set_axisbelow(True)
    # the log minor ticks previously rendered as a solid black band along the
    # bottom spine; keep a plain outward major tick only
    ax_c.minorticks_off()
    ax_c.set_xticks([200, 1000, 4000])
    ax_c.set_xlim(xlo, xhi * 2.35)
    ax_c.tick_params(axis="x", which="both", direction="out", length=3.2, width=0.7)
    ax_c.tick_params(axis="y", which="both", length=0)

    # ---------------- panel letters, drop the top and right spines ---------- #
    for ax, letter, ylab in zip([ax_b, ax_c], ["(b)", "(c)"],
                                ["Node age (Ma)", "Time factor"]):
        ax.text(PANEL_XY[0], PANEL_XY[1], letter, transform=ax.transAxes,
                fontsize=FS_PANEL, fontweight="bold", ha="left", va="bottom")
        ax.text(1.0, PANEL_XY[1], ylab, transform=ax.transAxes, fontsize=FS_AXIS,
                ha="right", va="bottom")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    a_axes[0].text(PANEL_XY[0] - 0.30, PANEL_XY[1], "(a)", transform=a_axes[0].transAxes,
                   fontsize=FS_PANEL, fontweight="bold", ha="left", va="bottom")


    for ax in (ax_b, ax_c):
        ax.yaxis.set_label_coords(-0.19, 0.5)

    # ---------------- equal-width self-check -------------------------------- #
    _fw = fig.get_size_inches()[0] * 25.4
    _widths = [ax.get_position().width * _fw for ax in (a_axes[0], ax_b, ax_c)]
    print("  panel widths (mm)      " + "  ".join("%.2f" % v for v in _widths))
    assert max(_widths) - min(_widths) < 0.05, (
        "the three columns must be the same width: %s" % _widths)
    print("  [width check] PASS - all three columns are the same width")

    # ---------------- bbox self-check before writing ------------------------ #
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    w_px, h_px = fig.canvas.get_width_height()
    offenders = check_all_bboxes_inside(fig, renderer, w_px, h_px)

    fig.savefig(os.path.join(out_dir, STEM + ".png"), dpi=1200, facecolor="white")
    fig.savefig(os.path.join(out_dir, STEM + ".pdf"), facecolor="white")
    plt.close(fig)
    _set_pdf_metadata(os.path.join(out_dir, STEM + ".pdf"),
                      "OpenRelTime Figure S3 - effective bounds, sensitivity and example-file conflict")

    print("wrote %s.png/.pdf from %s" % (STEM, DATA_JSON))
    print("  panel (a) x-spans (caption F4-4): " + ", ".join(
        "%s %.2f Ma" % (lab, sp) for lab, sp in spans))
    for row in rows:
        print("  %-28s bounds %.2f-%.2f  point %.2f  clipped=%s" % (
            row["node"], row["eff_min"], row["eff_max"], row["point"],
            row.get("on_lower_bound") or row.get("on_upper_bound")))
    print("  implied global factors: %s (spread %.1fx)" % (
        {k: round(v, 1) for k, v in implied.items()}, spread))
    print("  sweep: %d nodes over crown ages %s (all through origin)" % (
        len(nodes), xsc))
    print("[bbox check] page canvas = %.0f x %.0f px @ fig dpi %g" % (w_px, h_px, fig.dpi))
    if offenders:
        print("[bbox check] FAIL - %d artist(s) exceed the page rectangle:" % len(offenders))
        for o in offenders:
            print("    ", o)
        raise SystemExit(1)
    print("[bbox check] PASS - all text and artwork bboxes lie inside the page rectangle")
    overlaps = check_text_overlaps(fig, renderer)
    if overlaps:
        print("[overlap check] FAIL - %d collision(s):" % len(overlaps))
        for kind, a, b, ox, oy in overlaps[:20]:
            print("   ", kind, a, "<->", b, "overlap %.1f x %.1f px" % (ox, oy))
        raise SystemExit(2)
    print("[overlap check] PASS - no text overlaps other text or the plotted data")
    intrusions = check_column_territories(fig, renderer)
    if intrusions:
        print("[territory check] FAIL - %d label(s) reach into a neighbouring column:"
              % len(intrusions))
        for key, txt, bx0, bx1, (tx0, tx1) in intrusions[:20]:
            print("    col x0=%.3f %-30s spans %.1f..%.1f px, territory %.1f..%.1f px"
                  % (key, txt, bx0, bx1, tx0, tx1))
        raise SystemExit(2)
    print("[territory check] PASS - every label stays inside its own column")


def plain_linear(ax):
    fmt = FuncFormatter(lambda v, _: ("%g" % v))
    ax.xaxis.set_major_formatter(fmt)


def check_all_bboxes_inside(fig, renderer, w_px, h_px, tol=0.6):
    """Check that every text and artwork bbox lies inside the page rectangle."""
    offenders = []

    def add(artist, desc):
        try:
            bb = artist.get_window_extent(renderer=renderer)
        except TypeError:
            bb = artist.get_window_extent()
        x0, y0, x1, y1 = bb.x0, bb.y0, bb.x1, bb.y1
        if x0 < -tol or y0 < -tol or x1 > w_px + tol or y1 > h_px + tol:
            offenders.append((desc, round(float(x0), 1), round(float(y0), 1),
                              round(float(x1), 1), round(float(y1), 1)))

    for t in fig.texts:
        add(t, "figtext:" + repr(t.get_text()[:24]))
    for ax in fig.axes:
        for t in ax.texts:
            add(t, "text:" + repr(t.get_text()[:24].replace("\n", "/")))
        if ax.xaxis.label.get_text():
            add(ax.xaxis.label, "xlabel:" + repr(ax.xaxis.label.get_text()[:24].replace("\n", "/")))
        if ax.yaxis.label.get_text():
            add(ax.yaxis.label, "ylabel:" + repr(ax.yaxis.label.get_text()[:24].replace("\n", "/")))
        if ax.get_title():
            add(ax.title, "title:" + repr(ax.get_title()[:24]))
        for lab in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            if lab.get_text():
                add(lab, "tick:" + repr(lab.get_text()))
        # artwork whose display box runs off the page counts as an offender
        for artist in list(ax.lines) + list(ax.patches) + list(ax.collections):
            add(artist, "artist:" + type(artist).__name__)
    return offenders


def _set_pdf_metadata(pdf_path, title):
    """Stamp readable XMP/DocInfo metadata onto the exported PDF (optional; needs fitz)."""
    try:
        import fitz
    except ImportError:
        return False
    doc = fitz.open(pdf_path)
    md = doc.metadata or {}
    md.update({"title": title, "author": "", "creator": "OpenRelTime",
               "producer": "matplotlib",
               "keywords": "OpenRelTime; Bioinformatics; Application Note; "
                           "supplementary figure"})
    doc.set_metadata(md)
    tmp = pdf_path + ".tmp"
    doc.save(tmp); doc.close()
    os.replace(tmp, pdf_path)
    return True


def _cli():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--engine-root", default=None,
                    help="folder of the OpenRelTime repository (the one holding "
                         "openreltime/ and data/examples); only --regen-data "
                         "needs it; also read from $%s" % ENV_ENGINE)
    ap.add_argument("--golden-dir", default=None,
                    help="accepted for symmetry with the S2 script; Figure S3 "
                         "reads no R3F golden file, so this is ignored")
    ap.add_argument("--outdir", default=None,
                    help="where to write the .pdf/.png (default: this folder)")
    ap.add_argument("--regen-data", action="store_true",
                    help="recompute the data with the installed openreltime "
                         "(method='effective', n_effective=%d, seed=%d) before "
                         "drawing; needs the package and the repository.  The "
                         "refresh writes %s, NOT the checked-in %s, unless "
                         "--force is given" % (N_EFFECTIVE, SEED, REGEN_JSON, DATA_JSON))
    ap.add_argument("--regen-out", default=None,
                    help="where --regen-data should write (default: %s beside "
                         "this script); use e.g. /tmp/s3-data.json to compare "
                         "against the published record" % REGEN_JSON)
    ap.add_argument("--force", action="store_true",
                    help="let --regen-data overwrite the checked-in %s (the "
                         "published record) with this machine's numbers" % DATA_JSON)
    args = ap.parse_args()
    if args.golden_dir:
        print("note: --golden-dir is not used by Figure S3 (no R3F golden input)")
    if args.regen_data:  # pragma: no cover
        regen_data(args.engine_root, args.regen_out, args.force)
    main(args.outdir)


if __name__ == "__main__":
    _cli()
