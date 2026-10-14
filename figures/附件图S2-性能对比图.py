# -*- coding: utf-8 -*-
"""OpenRelTime Supplementary Figure S2 -- agreement with R3F and wall-clock time.

Reads every number from the data file checked in beside this script
(``附件图S2-性能对比图-data.json``) plus the two R3F golden CSVs that ship in the
OpenRelTime repository (``data/golden/r3f/example_RRF_rates.csv`` and
``example_RRF_times.csv``).  The render path never imports ``openreltime``, so
the figure reproduces on a machine that does not have the package installed.
One run writes a PDF exactly 178 mm wide and a 1200 dpi PNG, then measures the
delivered files and checks that every text and artwork bounding box lies inside
the page rectangle.

Three panels:
  (a) residuals of the agreement with R3F: x = the R3F value (log axis, trimmed
      to the span actually present in the data), y = the OpenRelTime residual
      relative to R3F in units of double-precision epsilon, with the +/-1e-13
      band, the +/-epsilon lines and the slope / median |delta| / max |delta| /
      bit-identical counts for rates and node times.
  (b) wall-clock time against tree size (log-log): median line with min-max
      whiskers (visible repeat dispersion, in place of invisible IQR whiskers).
      R3F is NOT timed here: the only R3F mark is an explicit upper bound drawn
      at the "less than one second for 1,000 tips" value reported in the R3F
      publication, labelled as such inside the panel.
  (c) dispersion between repeats: the coefficient of variation at each size as
      one equal-width bar per size on a categorical axis.

Run:
    python 附件图S2-性能对比图.py
    python 附件图S2-性能对比图.py --golden-dir /path/to/data/golden/r3f
    python 附件图S2-性能对比图.py --engine-root /path/to/OpenRelTime

Refresh the data (optional; needs the OpenRelTime repository).  The checked-in
-data.json IS the published record, so the refresh writes a COPY beside it
(``附件图S2-性能对比图-data.regen.json``) unless you name a file with
``--regen-out PATH`` or consciously pass ``--force``:

    python 附件图S2-性能对比图.py --regen-data --engine-root /path/to/OpenRelTime
    python 附件图S2-性能对比图.py --regen-data --regen-out /tmp/s2-data.json
    python 附件图S2-性能对比图.py --regen-data --force     # replaces -data.json

Read that guard as part of the data-availability statement: ``--regen-data``
re-measures the wall-clock benchmark on the machine it runs on, so it is the way
to check the agreement statistics of panel (a), never a way to verify the
published timings of panels (b) and (c), which belong to the machine recorded
in the file's ``environment`` block.  ``--force`` exists for a deliberate
re-measurement that is meant to become the new record.
"""

import argparse
import json
import os
import platform
import subprocess
import sys

import numpy as np
from matplotlib.transforms import Bbox
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
PROJECT_ROOT = os.path.dirname(HERE)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import LogLocator, NullFormatter  # noqa: E402

STEM = "附件图S2-性能对比图"
DATA_JSON = STEM + "-data.json"
# --regen-data writes here unless --force is given: -data.json is the published
# record of the timings in Table S5, and a refresh on another machine measures
# that machine, not the record (review item S2-i / F3-7).
REGEN_JSON = STEM + "-data.regen.json"
GOLDEN_FILES = ("example_RRF_rates.csv", "example_RRF_times.csv")
GOLDEN_SUBPATH = os.path.join("data", "golden", "r3f")
ENV_GOLDEN = "OPENRELTIME_GOLDEN_DIR"
ENV_ENGINE = "OPENRELTIME_ENGINE_ROOT"


def _looks_like_engine(path):
    return (os.path.isdir(os.path.join(path, "openreltime"))
            and os.path.isfile(os.path.join(path, GOLDEN_SUBPATH, GOLDEN_FILES[0])))


def _candidate_dirs(extra):
    """Folders worth probing for the repository, in priority order."""
    out = []
    for given in (extra, os.environ.get(ENV_ENGINE), PROJECT_ROOT, os.getcwd()):
        if given:
            out.append(os.path.abspath(os.path.expanduser(given)))
    siblings = PROJECT_ROOT
    try:
        out += [os.path.join(siblings, n) for n in sorted(os.listdir(siblings))]
    except OSError:
        pass
    return out


def locate_engine_root(explicit=None):
    """The OpenRelTime repository folder, found by content rather than by name.

    The repository's folder name is deliberately not assumed: as shipped to a
    journal the sibling folders may be renamed or absent, so every candidate is
    probed for ``openreltime/`` plus ``data/golden/r3f/``.
    """
    for path in _candidate_dirs(explicit):
        if _looks_like_engine(path):
            return path
    return None


def locate_golden_dir(explicit=None, engine_root=None):
    """The folder holding the two R3F golden CSVs used by panel (a)."""
    for given in (explicit, os.environ.get(ENV_GOLDEN)):
        if given:
            path = os.path.abspath(os.path.expanduser(given))
            if all(os.path.isfile(os.path.join(path, f)) for f in GOLDEN_FILES):
                return path
    roots = [engine_root] if engine_root else []
    roots += [p for p in _candidate_dirs(None) if _looks_like_engine(p)]
    for root in roots:
        if not root:
            continue
        path = os.path.join(root, GOLDEN_SUBPATH)
        if all(os.path.isfile(os.path.join(path, f)) for f in GOLDEN_FILES):
            return path
    return None


def resolve_paths(golden_dir=None, engine_root=None):
    """Return (golden_dir, engine_root) or stop with a message naming the flag."""
    root = locate_engine_root(engine_root)
    golden = locate_golden_dir(golden_dir, root)
    if golden is None:
        raise SystemExit(
            "error: the R3F golden CSVs (%s) were not found.\n"
            "  Panel (a) plots OpenRelTime against those files, so they are\n"
            "  needed to draw the figure.  Point this script at them with\n"
            "  --golden-dir PATH (a folder containing %s), or at the OpenRelTime\n"
            "  repository with --engine-root PATH, or set %s / %s."
            % (", ".join(GOLDEN_FILES), "/".join(GOLDEN_FILES), ENV_GOLDEN, ENV_ENGINE))
    return golden, root

# --- shared style constants, identical to the S3 script -------------------- #
MM = 1.0 / 25.4
COLUMN_MM = 178.0                      # double-column text block: the figure is
                                        # drawn at final size, never rescaled
FS_PANEL = 9.0                         # panel letters
FS_AXIS = 8.5                          # axis labels
FS_TICK = 8.0                          # tick labels
FS_SERIES = 7.0                        # series and in-panel labels
FS_NOTE = 7.5                          # in-panel notes: the engine's own floor
                                       # (make_figures.py MIN_PT = 7.5) is what
                                       # the statistics box and the R3F bound
                                       # note were undercutting at 6.0/6.5 pt
PANEL_XY = (-0.02, 1.02)               # panel-letter position, same in both figures
# Okabe-Ito / colour-blind-safe palette
C_DATA = "#0072B2"                     # OpenRelTime measurements (main blue, as S3)
C_DATA2 = "#E69F00"                    # second measured series (orange)
C_REF = "#D55E00"                      # R3F reference / literature value (vermilion)
C_BAND = "#56B4E9"                     # residual reference band (pale blue, shading)
C_EPS = "#009E73"                      # epsilon reference line (green)
C_GREY = "#555555"
MACHINE_EPS = 2.220446049250313e-16
# F3-1 / P1-20: the panel must not read as "sub-epsilon agreement".  The two
# golden CSVs are printed at 15 significant digits (%.15g), so their own
# quantisation is 2.25-22.5 machine epsilons depending on the leading digit.
REFERENCE_NOTE = ("agreement is against the\n"
                  "15-s.d. printed reference\n"
                  "(its floor: 2.25-22.5 ε)")

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
    "font.size": FS_TICK, "axes.labelsize": FS_AXIS,
    "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
    "axes.linewidth": 0.8, "axes.edgecolor": "#333333",
    "savefig.dpi": 1200, "figure.dpi": 100,
})


def plain_log_ticks(axis):
    """Label a log axis with plain decimal numbers (0.001 / 0.01 / 1 ...).

    This removes the 5.6 pt 10^n superscripts.  ``axis`` must be an Axes
    xaxis/yaxis instance and this must be called *after* the view limits are
    set, so that only the decades inside the range are labelled and no out of
    range decade label lands on the page edge.
    """
    import math
    from matplotlib.ticker import FuncFormatter
    vmin, vmax = axis.get_view_interval()
    lo = int(math.floor(math.log10(max(vmin, 1e-12))))
    hi = int(math.ceil(math.log10(vmax)))
    ticks = [10.0 ** k for k in range(lo, hi + 1) if vmin <= 10.0 ** k <= vmax]
    axis.set_ticks(ticks)
    axis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    # S2-e: this installed 90 minors per decade, which at 178 mm rendered as a
    # solid black band along both spines of panels (a) and (b).  Two minors per
    # decade (0.2 and 0.5) keep the log scale readable; the sibling S3 script
    # documents the same fix at its own panel (c).
    axis.set_minor_locator(LogLocator(base=10.0, subs=[0.2, 0.5], numticks=100))
    axis.set_minor_formatter(NullFormatter())


# --------------------------------------------------------------------------- #
# pandas/numpy helpers: rebuild the clade map from the R3F golden CSVs without
# importing openreltime
# --------------------------------------------------------------------------- #
def clade_map(df):
    clades, pending, rows = {}, [], {}
    for record in df.to_dict("records"):
        nid = int(record["NodeId"])
        rows[nid] = record
        if str(record["NodeLabel"]) != "-":
            clades[nid] = frozenset([record["NodeLabel"]])
        else:
            pending.append(nid)
    while pending:
        progressed = False
        for nid in list(pending):
            row = rows[nid]
            try:
                d1, d2 = int(row["Des1"]), int(row["Des2"])
            except (TypeError, ValueError):
                continue
            if d1 in clades and d2 in clades:
                clades[nid] = clades[d1] | clades[d2]
                pending.remove(nid)
                progressed = True
        if not progressed:
            break
    return clades


def golden_by_clade(df, column, tips_as_zero=False):
    clades = clade_map(df)
    values = {}
    for record in df.to_dict("records"):
        nid = int(record["NodeId"])
        if str(record["NodeLabel"]) != "-":
            values[frozenset([record["NodeLabel"]])] = (
                0.0 if tips_as_zero else float(record[column]))
        else:
            values[clades[nid]] = float(record[column])
    return values


def _canon(clade):
    return "|".join(sorted(clade))


def agreement_points(payload, golden_dir):
    """x = the R3F value (from the golden CSVs), y = the OpenRelTime value held in
    the per-node arrays of the checked-in -data.json."""
    arrays = payload["agreement_arrays"]
    orr, ort = arrays["openreltime_rate"], arrays["openreltime_time"]
    g_rates = golden_by_clade(pd.read_csv(os.path.join(golden_dir, "example_RRF_rates.csv")), "Rate")
    g_times = golden_by_clade(pd.read_csv(os.path.join(golden_dir, "example_RRF_times.csv")),
                              "Time", tips_as_zero=True)
    rx, ry = [], []
    for clade, v in g_rates.items():
        k = _canon(clade)
        if v > 0 and k in orr:
            rx.append(v); ry.append(orr[k])
    tx, ty = [], []
    for clade, v in g_times.items():
        k = _canon(clade)
        if v > 0 and k in ort:
            tx.append(v); ty.append(ort[k])
    return (np.array(rx), np.array(ry)), (np.array(tx), np.array(ty))


# --------------------------------------------------------------------------- #
# Only --regen-data imports openreltime and re-runs the benchmark, overwriting
# -data.json.  The default render path never enters this function, which is what
# keeps the figure reproducible without the package installed.
# --------------------------------------------------------------------------- #
def _sysctl(key: str):
    """Read one ``sysctl`` string/int on macOS; ``None`` if unavailable."""
    try:
        out = subprocess.run(["sysctl", "-n", key], capture_output=True,
                             text=True, timeout=10, check=False).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None
    if not out:
        return None
    return int(out) if out.isdigit() else out


def _machine_label() -> str:
    """Human-usable machine identity, e.g. ``Apple M5 (Mac16,1, arm64)``."""
    model = _sysctl("hw.model")
    arch = platform.machine()
    bits = [x for x in (_sysctl("machdep.cpu.brand_string"), model, arch) if x]
    return " / ".join(str(b) for b in bits) or platform.platform()


def regen_data(engine_root=None, golden_dir=None, out_json=None,
               force=False):  # pragma: no cover - needs openreltime
    import statistics
    import time
    # Resolve the repository FIRST and put it on sys.path: --engine-root points
    # at a source checkout, not an installed distribution, so importing
    # openreltime before the lookup failed even when the flag was correct.
    root = locate_engine_root(engine_root)
    if root is None:
        raise SystemExit(
            "error: --regen-data needs the OpenRelTime repository (it reads\n"
            "  data/examples/example.nwk and data/golden/r3f/).  Pass\n"
            "  --engine-root PATH, or set %s." % ENV_ENGINE)
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    import openreltime
    from openreltime import PhyloNode, read_tree, rrf_rates, rrf_times
    out_json = out_json or (DATA_JSON if force else REGEN_JSON)
    data_dir = os.path.join(root, "data")
    ex_tree = os.path.join(data_dir, "examples", "example.nwk")
    OUTGROUP = ["Ornithorhynchus_anatinus", "Zaglossus_bruijni", "Tachyglossus_aculeatus"]
    SEED, WARMUP = 20260906, 3
    SIZES = [100, 200, 500, 1000, 2000, 5000, 10000]
    REPEATS = {100: 51, 200: 51, 500: 31, 1000: 31, 2000: 21, 5000: 15, 10000: 9}

    def random_tree(n, rng):
        sub = [PhyloNode(label="T%d" % i, blen=float(rng.uniform(0.01, 0.5))) for i in range(n)]
        while len(sub) > 1:
            idx = rng.permutation(len(sub)); merged = []
            for k in range(0, len(idx) - 1, 2):
                p = PhyloNode(blen=float(rng.uniform(0.01, 0.5)))
                p.add_child(sub[idx[k]], blen=sub[idx[k]].blen)
                p.add_child(sub[idx[k + 1]], blen=sub[idx[k + 1]].blen)
                merged.append(p)
            if len(idx) % 2 == 1:
                merged.append(sub[idx[-1]])
            sub = merged
        return sub[0]

    tree = read_tree(ex_tree, outgroup=OUTGROUP)
    my_rates, my_times = rrf_rates(tree), rrf_times(tree)
    cl_of = {n.node_id: n.clade for n in tree.walk()}
    m_rates = {cl_of[i]: v for i, v in my_rates.rates.items()}
    m_times = {cl_of[i]: v for i, v in my_times.times.items()}
    golden = locate_golden_dir(golden_dir, root)
    g_rates = golden_by_clade(pd.read_csv(os.path.join(golden, "example_RRF_rates.csv")), "Rate")
    g_times = golden_by_clade(pd.read_csv(os.path.join(golden, "example_RRF_times.csv")),
                              "Time", tips_as_zero=True)
    rc = sorted(set(g_rates) & set(m_rates), key=lambda c: sorted(c))
    tc = sorted(set(g_times) & set(m_times), key=lambda c: sorted(c))
    bench = []
    for size in SIZES:
        warm = random_tree(size, np.random.default_rng(SEED))
        for _ in range(WARMUP):
            rrf_times(warm)
        ts = []
        for k in range(REPEATS[size]):
            cand = random_tree(size, np.random.default_rng(SEED + 1000 * size + k))
            t0 = time.perf_counter(); rrf_times(cand); ts.append(time.perf_counter() - t0)
        bench.append({"tips": size, "reps": len(ts), "median_s": statistics.median(ts),
                      "min_s": min(ts), "max_s": max(ts),
                      "q1_s": float(np.percentile(ts, 25)), "q3_s": float(np.percentile(ts, 75)),
                      "cv_pct": 100.0 * statistics.pstdev(ts) / statistics.mean(ts)})
    with open(DATA_JSON, "r", encoding="utf-8") as fh:
        payload = json.load(fh)
    payload["benchmark"] = bench
    payload["agreement_arrays"] = {
        "openreltime_rate": {_canon(c): float(m_rates[c]) for c in rc},
        "openreltime_time": {_canon(c): float(m_times[c]) for c in tc},
        "n_rate": len(rc), "n_time": len(tc)}
    # Record the measurement environment in the file itself.  The published
    # timings were produced before this block existed, and the manuscript now
    # discloses that interpreter from these keys rather than from prose.
    payload["openreltime_version"] = getattr(openreltime, "__version__", "unknown")
    payload["python"] = platform.python_version()
    payload["numpy"] = np.__version__
    payload["pandas"] = pd.__version__
    payload["matplotlib"] = matplotlib.__version__
    payload["environment"] = {
        "machine": _machine_label(),
        "os_version": "macOS %s (%s)" % (platform.mac_ver()[0], platform.mac_ver()[2]),
        "python": platform.python_version(),
        "clock": "time.perf_counter",
        "cpu_model": _sysctl("machdep.cpu.brand_string"),
        "physical_memory_bytes": _sysctl("hw.memsize"),
    }
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
    if os.path.abspath(out_json) == os.path.abspath(DATA_JSON):
        print("refreshed (forced, published record replaced)", out_json)
    else:
        print("refreshed a COPY at", out_json, "-- the published record",
              DATA_JSON, "is untouched; pass --force to replace it")


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
            offenders.append((desc, round(x0, 1), round(y0, 1), round(x1, 1), round(y1, 1)))

    for t in fig.texts:
        add(t, "figtext:" + repr(t.get_text()[:24]))
    for ax in fig.axes:
        for t in ax.texts:
            add(t, "text:" + repr(t.get_text()[:24].replace("\n", "/")))
        add(ax.xaxis.label, "xlabel:" + repr(ax.xaxis.label.get_text()[:24].replace("\n", "/")))
        add(ax.yaxis.label, "ylabel:" + repr(ax.yaxis.label.get_text()[:24].replace("\n", "/")))
        if ax.get_title():
            add(ax.title, "title:" + repr(ax.get_title()[:24]))
        for lab in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            if lab.get_text():
                add(lab, "tick:" + repr(lab.get_text()))
        # artwork: lines, scatter points, bars and whisker boxes must stay on the page too
        for artist in list(ax.lines) + list(ax.patches) + list(ax.collections):
            add(artist, "artist:" + type(artist).__name__)
    return offenders


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


def main(golden_dir=None, engine_root=None, out_dir=None):
    out_dir = HERE if out_dir is None else out_dir
    golden = locate_golden_dir(golden_dir, locate_engine_root(engine_root))
    with open(DATA_JSON, "r", encoding="utf-8") as fh:  # lives beside the script
        payload = json.load(fh)
    a = payload["agreement_with_r3f"]
    bench = payload["benchmark"]
    (rx, ry), (tx, ty) = agreement_points(payload, golden)

    # F3-2 / S2-h: the stored scalars of agreement_with_r3f are quoted in the
    # in-panel box next to the points they describe, and the shipped file's
    # time_slope (0.99999999999999944) is not the slope of the shipped arrays
    # (numpy.polyfit gives 0.99999999999999978).  Derive every quoted scalar
    # from the arrays that are actually plotted, so the box can never disagree
    # with them.  The stored arrays are read-only here, and every recomputed
    # value agrees with the stored one to the printed digit (the slopes print
    # 1.000000 either way), so the published numbers do not move.
    def _stats(x, y):
        rel = np.abs((y - x) / x)
        return {"slope": float(np.polyfit(x, y, 1)[0]),
                "med": float(np.median(rel)), "max": float(np.max(rel)),
                "bit": int(np.sum(x == y)), "n": int(len(x))}
    rstat, tstat = _stats(rx, ry), _stats(tx, ty)
    for key, got, want in (("rate_slope", rstat["slope"], a["rate_slope"]),
                           ("rate_median_rel", rstat["med"], a["rate_median_rel"]),
                           ("rate_max_rel", rstat["max"], a["rate_max_rel"]),
                           ("rate_bit_identical", rstat["bit"], a["rate_bit_identical"]),
                           ("n_rate_entries", rstat["n"], a["n_rate_entries"]),
                           ("time_slope", tstat["slope"], a["time_slope"]),
                           ("time_median_rel", tstat["med"], a["time_median_rel"]),
                           ("time_max_rel", tstat["max"], a["time_max_rel"]),
                           ("time_bit_identical", tstat["bit"], a["time_bit_identical"]),
                           ("n_time_positive", tstat["n"], a["n_time_positive"])):
        if got != want:
            print("  [scalars] %-18s stored %r -> recomputed %r (from the plotted arrays)"
                  % (key, want, got))

    fig, axes = plt.subplots(1, 3, figsize=(COLUMN_MM * MM, COLUMN_MM * MM * 0.50))
    ax_a, ax_b, ax_c = axes

    # ---------------- (a) residuals vs the printed R3F reference ------------ #
    eps = MACHINE_EPS
    rres = (ry - rx) / rx / eps          # relative residual in units of eps
    tres = (ty - tx) / tx / eps
    band = 1e-13 / eps                    # how many eps make +/-1e-13
    ax_a.axhspan(-band, band, color=C_BAND, alpha=0.22, lw=0, zorder=0)
    for s in (1, -1):
        ax_a.axhline(s, color=C_EPS, lw=0.9, ls="--", zorder=1)
    ax_a.axhline(0.0, color=C_GREY, lw=0.7, zorder=1)
    ax_a.scatter(rx, rres, s=6, color=C_DATA, marker="o", linewidths=0,
                 alpha=0.85, zorder=3, label="rates (n=%d)" % a["n_rate_entries"])
    ax_a.scatter(tx, tres, s=6, color=C_DATA2, marker="s", linewidths=0,
                 alpha=0.85, zorder=3, label="node times (n=%d)" % a["n_time_positive"])
    ax_a.set_xscale("log")
    ax_a.set_xlim(min(rx.min(), tx.min()) * 0.8, max(rx.max(), tx.max()) * 1.25)
    # the axis range leaves 15% of margin over the +/-1e-13 band; the statistics
    # box is inside the panel and quotes the caption word for word
    ax_a.set_ylim(-band * 1.15, band * 1.15)
    ax_a.set_xlabel("R3F value (rate or node time)")
    ax_a.set_ylabel("Rel. residual / \u03b5")
    plain_log_ticks(ax_a.xaxis)
    # The statistics box repeats, verbatim, the slope / median |delta| / max
    # |delta| / bit-identical counts given in the caption, but every number is
    # now derived above from the very arrays that are plotted.  It sits in the
    # empty top-right of the panel (the top-left is taken by the largest
    # residual) and is re-flowed into short lines so that the 7.5 pt type (the
    # engine's own floor; it was 6.0 pt) cannot reach the neighbouring column.
    # F3-1 / P1-20: the box is quoted word for word in the caption, so it stays
    # eight lines; the reference-floor qualifier that the caption now carries
    # goes into its own note in the empty band below the points (every residual
    # lies between 0.48 and 0.66 of the panel height, so the strip under 0.45 is
    # free).  Both were 6.0 pt, under the engine's own 7.5 pt floor (F3-6/F6-5).
    stats_note = (
        "rates: slope = %.6f\n"
        "  median |Δ| = %.2f ε\n"
        "  max |Δ| = %.1f ε\n"
        "  %d of %d bit-identical\n"
        "times: slope = %.6f\n"
        "  median |Δ| = %.2f ε\n"
        "  max |Δ| = %.1f ε\n"
        "  %d of %d bit-identical"
        % (rstat["slope"], rstat["med"] / eps, rstat["max"] / eps,
           rstat["bit"], rstat["n"],
           tstat["slope"], tstat["med"] / eps, tstat["max"] / eps,
           tstat["bit"], tstat["n"])
    )
    ax_a.text(0.975, 0.975, stats_note, transform=ax_a.transAxes, fontsize=FS_NOTE,
              va="top", ha="right", zorder=5, linespacing=1.15,
              bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#bbbbbb",
                        lw=0.5, alpha=0.92))
    ax_a.text(0.02, 0.42, REFERENCE_NOTE, transform=ax_a.transAxes, fontsize=FS_NOTE,
              va="top", ha="left", zorder=5, linespacing=1.15, color="#333333")
    leg = ax_a.legend(loc="lower right", fontsize=FS_SERIES, framealpha=0.9,
                      borderpad=0.4, labelspacing=0.45)
    for h, lab in zip(leg.legend_handles, leg.get_texts()):
        lab.set_color(h.get_facecolor())

    # ---------------- (b) wall-clock time against tree size ----------------- #
    tips = np.array([r["tips"] for r in bench], float)
    med = np.array([r["median_s"] for r in bench])
    lo = med - np.array([r["min_s"] for r in bench])
    hi = np.array([r["max_s"] for r in bench]) - med
    # S2-d / F3-4: the line is the median of the repeats and the whiskers are
    # their min and max, which the caption said but the panel did not.
    ax_b.errorbar(tips, med, yerr=[lo, hi], fmt="o-", color=C_DATA, ms=4.2, lw=1.2,
                  elinewidth=1.0, capsize=2.8, capthick=1.0, zorder=3,
                  label="median (n repeats)")
    ax_b.legend(loc="lower right", fontsize=FS_SERIES, framealpha=0.9,
                borderpad=0.4, labelspacing=0.45)
    ax_b.set_xscale("log"); ax_b.set_yscale("log")
    ax_b.set_xlim(80, 13000); ax_b.set_ylim(2.2e-4, 2.6)   # headroom for the bound
    ax_b.set_xlabel("Number of tips")
    ax_b.set_ylabel("")          # moved above the axes: no label in the gutter
    plain_log_ticks(ax_b.xaxis)
    plain_log_ticks(ax_b.yaxis)
    # R3F is not timed here.  Its publication reports that the RRF functions
    # "took less than one second for a phylogeny with 1000 tips", averaged over
    # rrf_rates_times, corrtest and ddbd on two simulated trees -- an upper
    # bound on a different quantity from the rrf_times-only benchmark above, so
    # it is drawn as a bound (dashed line, hollow square, downward arrow) and
    # labelled as unmeasured, never as a point estimate.
    ax_b.hlines(1.0, 350, 3000, colors=C_REF, lw=1.0, ls=(0, (3.0, 2.0)), zorder=4)
    ax_b.plot([1000], [1.0], marker="s", ms=8, mfc="none", mec=C_REF, mew=1.6,
              ls="none", zorder=5)   # hollow vermilion square: the reported bound
    ax_b.annotate("", xy=(1000, 0.55), xytext=(1000, 0.90), zorder=5,
                  arrowprops=dict(arrowstyle="-|>", color=C_REF, lw=1.0,
                                  shrinkA=0, shrinkB=0, mutation_scale=9.0))
    ax_b.text(1200, 0.90, "R3F: reported\n<1 s at 1,000 tips\n(not measured here)",
              fontsize=FS_NOTE, color=C_REF, ha="left", va="top", linespacing=1.3,
              zorder=5)

    # ---------------- (c) dispersion between repeats ------------------------ #
    cv = np.array([r["cv_pct"] for r in bench])
    slots = np.arange(len(cv))
    ax_c.bar(slots, cv, width=0.66, color=C_DATA, alpha=0.9,
             edgecolor="#004d78", lw=0.5, zorder=2)
    ax_c.set_xlim(-0.62, len(cv) - 0.38)
    ax_c.set_xticks(slots)
    ax_c.set_xticklabels(["%d" % int(v) for v in tips], rotation=90, ha="center",
                         va="top", fontsize=FS_TICK)
    ax_c.set_ylim(0, cv.max() * 1.12)
    ax_c.grid(axis="y", lw=0.4, color="#ededed", zorder=0)
    ax_c.set_axisbelow(True)
    ax_c.set_xlabel("Number of tips")
    ax_c.set_ylabel("")


    # ---------------- panel letters, drop the top and right spines ---------- #
    # S2-d: panel (b)'s title names the statistic, not just the quantity; it is
    # wrapped in two lines so the wider title cannot reach into panel (a).
    for ax, letter, ylab in zip((ax_a, ax_b, ax_c), ["(a)", "(b)", "(c)"],
                                ["Rel. residual / \u03b5",
                                 "Wall-clock time (s),\nmedian with min\u2013max",
                                 "CV of repeat times (%)"]):
        ax.text(PANEL_XY[0], PANEL_XY[1], letter, transform=ax.transAxes,
                fontsize=FS_PANEL, fontweight="bold", ha="left", va="bottom")
        if ax is not ax_a:
            # the y title sits on the panel's own top edge, right aligned, so it
            # can never reach into the neighbouring column through the gutter
            ax.text(1.0, PANEL_XY[1], ylab, transform=ax.transAxes, fontsize=FS_AXIS,
                    ha="right", va="bottom")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.subplots_adjust(left=0.090, right=0.986, bottom=0.205, top=0.882, wspace=0.52)
    for ax in (ax_a, ax_b, ax_c):   # fixed y-label position: the gutter is designed
        ax.yaxis.set_label_coords(-0.26, 0.5)
        # S2-e: 1.8 pt minors over 4 decades read as a bar; 1.2 pt with two
        # minors per decade (see plain_log_ticks) does not
        ax.tick_params(axis="both", which="minor", length=1.2, width=0.5)


    # ---------------- bbox self-check before writing ------------------------ #
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    w_px, h_px = fig.canvas.get_width_height()
    offenders = check_all_bboxes_inside(fig, renderer, w_px, h_px)

    png = os.path.join(out_dir, STEM + ".png")
    pdf = os.path.join(out_dir, STEM + ".pdf")
    fig.savefig(png, dpi=1200, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)
    _set_pdf_metadata(pdf,
                      "OpenRelTime Figure S2 - agreement with R3F and runtime")

    print("wrote %s.png/.pdf from %s (+ golden CSVs)" % (STEM, DATA_JSON))
    print("  agreement: rate max %.2e  time max %.2e  (%d+%d bit-identical of %d+%d)" % (
        a["rate_max_rel"], a["time_max_rel"], a["rate_bit_identical"], a["time_bit_identical"],
        a["n_rate_entries"], a["n_time_positive"]))
    for r in bench:
        print("  %6d tips  median %.5f s  min-max %.5f-%.5f  CV %.1f%%  (n=%d)" % (
            r["tips"], r["median_s"], r["min_s"], r["max_s"], r["cv_pct"], r["reps"]))
    print("[bbox check] page canvas = %.0f x %.0f px @ fig dpi %g" % (w_px, h_px, fig.dpi))
    if offenders:
        print("[bbox check] FAIL - %d artist(s) exceed the page rectangle:" % len(offenders))
        for o in offenders:
            print("    ", o)
        sys.exit(1)
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
                         "openreltime/ and data/golden/r3f); also read from $%s"
                         % ENV_ENGINE)
    ap.add_argument("--golden-dir", default=None,
                    help="folder holding %s; also read from $%s"
                         % (" and ".join(GOLDEN_FILES), ENV_GOLDEN))
    ap.add_argument("--outdir", default=None,
                    help="where to write the .pdf/.png (default: this folder)")
    ap.add_argument("--regen-data", action="store_true",
                    help="re-run the benchmark and the agreement comparison with "
                         "the installed openreltime before drawing (needs the "
                         "package; not needed to draw the figure).  Writes %s, "
                         "NOT the checked-in %s, unless --force is given"
                         % (REGEN_JSON, DATA_JSON))
    ap.add_argument("--regen-out", default=None,
                    help="where --regen-data should write (default: %s beside "
                         "this script)" % REGEN_JSON)
    ap.add_argument("--force", action="store_true",
                    help="let --regen-data overwrite the checked-in %s with this "
                         "machine's timings (the published record; see the "
                         "module docstring)" % DATA_JSON)
    args = ap.parse_args()
    if args.regen_data:  # pragma: no cover
        regen_data(args.engine_root, args.golden_dir, args.regen_out, args.force)
    main(args.golden_dir, args.engine_root, args.outdir)


if __name__ == "__main__":
    _cli()
