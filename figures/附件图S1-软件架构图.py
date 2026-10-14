# -*- coding: utf-8 -*-
"""OpenRelTime Supplementary Figure S1 -- software architecture.

Layout and wording only; every drawing primitive, the printing constraints, the
PDF/PNG export, the delivery checks and the editable ``.drawio`` mirror come
from ``make_figures.py``, which also draws Figure 1.  S1 and Figure 1 therefore
share one type size, one line weight, one colour ramp and one shape grammar.

The figure shows five lanes: the user interfaces, the three core engines (rrf,
with its public rates and times wrappers, calibrate, ci), the input/output
layer, the optional extension modules and the foundation.  Every module of the
released distribution is named on a card, and each arrow is a real dependency
or data flow, not decoration: a dashed arrow leaves the module that imports and
points at the module that is imported, and the single solid arrow follows the
value that is passed at run time.

Run:
    python 附件图S1-软件架构图.py
    python 附件图S1-软件架构图.py --figure-engine /path/to/folder/with/make_figures.py

``make_figures.py`` is looked up in ``--figure-engine``, in
$OPENRELTIME_FIGURE_ENGINE, in this folder, and in any sibling folder that
contains it; the folder's name is never assumed.
"""

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIGURE_ENGINE = "make_figures.py"
ENV_VAR = "OPENRELTIME_FIGURE_ENGINE"


def locate_figure_engine(explicit=None):
    """Find the folder holding ``make_figures.py`` without naming any folder.

    Order: ``--figure-engine``, ``$OPENRELTIME_FIGURE_ENGINE``, this folder, the
    sibling folders of this folder, the working directory.  Returns ``None`` when
    nothing matches, so the caller can raise a message that says what to set.
    """
    def candidate(path):
        path = os.path.abspath(os.path.expanduser(path))
        return path if os.path.isfile(os.path.join(path, FIGURE_ENGINE)) else None

    for given in (explicit, os.environ.get(ENV_VAR), HERE, os.getcwd()):
        if given:
            found = candidate(given)
            if found:
                return found
    siblings = os.path.abspath(os.path.join(HERE, os.pardir))
    try:
        names = sorted(os.listdir(siblings))
    except OSError:
        names = []
    for name in names:
        found = candidate(os.path.join(siblings, name))
        if found:
            return found
    return None


def _import_engine(explicit=None):
    folder = locate_figure_engine(explicit)
    if folder is None:
        raise SystemExit(
            "error: cannot find %s.\n"
            "  Figure S1 shares its drawing engine with Figure 1.  Point this\n"
            "  script at the folder that holds it with --figure-engine PATH, or\n"
            "  set %s to that folder, or run both figures from the same\n"
            "  folder." % (FIGURE_ENGINE, ENV_VAR))
    sys.path.insert(0, folder)
    return __import__("make_figures")


def build(out_dir=None, engine_folder=None):
    """Render Figure S1 and return its measured delivery report."""
    M = _import_engine(engine_folder)
    from make_figures import (Card, PAGE_W, FRAME_INSET, PAD, LINE, FONT_BODY,
                              GREY_TXT, RAMP, RAMP_EDGE, W_STAGE, W_DEP, LINE_MIN,
                              text_mm, Sheet, verify, write_drawio, DASH_STYLE)

    # ------------------------------------------------------------------ #
    # grid: one 3-column track shared by every lane, equal box heights
    # ------------------------------------------------------------------ #
    FX0, FX1 = FRAME_INSET, PAGE_W - FRAME_INSET
    LBL_W = 23.0                # lane label cell: keeps the arrow channels clear
    CELL_PAD = 1.6
    BOX_GAP = 4.5
    AREA_X0 = FX0 + 1.6 + LBL_W + 2.6
    AREA_X1 = FX1 - 1.6
    COL_W = (AREA_X1 - AREA_X0 - 2 * BOX_GAP) / 3.0
    COL_X = [AREA_X0 + i * (COL_W + BOX_GAP) for i in range(3)]
    N_LINES = 3                 # name + two detail lines, every box the same
    BOX_H = N_LINES * LINE + 2 * PAD
    LANE_PAD = 2.0
    LANE_GAP = 4.0
    BUS_GAP = 7.4
    ROW_GAP = 1.7

    MOD = dict(group="main", bold_first=True)
    OPT = dict(group="optional", bold_first=True)

    # The engine's own optional extras are matplotlib (``[plot]``) and the
    # development tools (``[dev]``); PySide6 is not one of them -- it is a
    # required dependency of the separate OpenRelTime-Studio distribution.
    LANES = [
        ("User interfaces", [[
            Card(["Python API", "import openreltime", ""], **MOD),
            Card(["Command line", "openreltime / ort", ""], **MOD),
            Card(["Desktop GUI (Studio)", "separate distribution;",
                  "PySide6 required there"], **OPT),
        ]], 0),
        ("Core engines", [[
            # `rrf.py` holds the class; `rates.py` and `times.py` hold the
            # `rrf_rates` / `rrf_times` / `rrf_rates_times` functions that the
            # API (`__init__.py:63,79`) and the CLI (`cli.py:165,195,226`) call,
            # so all three modules of the RRF layer are named on this card.
            Card(["rrf · rates · times", "RrfEngine + wrappers,",
                  "postorder + preorder"], **MOD),
            Card(["calibrate", "factor f, violations,", "effective bounds"], **MOD),
            Card(["ci", "delta-method confidence", "intervals (ref. 4)", ""], **MOD),
        ]], 1),
        ("Input / output", [[
            Card(["treeio", "Newick / NEXUS I/O,", "[&rate=...] preserved", ""], **MOD),
            Card(["_external", "validated argument list,", "subprocesses (no shell)",
                  ""], **MOD),
            Card(["report", "CSV / NEXUS / JSON writers,", "shared by API and CLI"],
                 **MOD),
        ]], 2),
        ("Optional extension modules", [
            [Card(["table · taxonomy", "tree2table records,", "monophyly check"], **OPT),
             Card(["megacc bridge", "MEGA-CC interoperability", "(external program)",
                   ""], **OPT),
             Card(["viz", "timetree and interval plots", "(matplotlib, optional)", ""],
                  **OPT)],
            [Card(["ddbd", "birth-death tree prior", ""], **OPT),
             Card(["corrtest", "rate-autocorrelation test,", "fixed coefficients"], **OPT),
             Card(["bootstrap", "bag-of-little-bootstraps,",
                   "IQ-TREE via _external"], **OPT)],
        ], 3),
        ("Foundation", [
            [Card(["tree", "PhyloNode binary tree,", "postorder / preorder"], **MOD),
             Card(["_constants", "epsilons, thresholds,", "CorrTest coefficients"], **MOD),
             Card(["third-party stack", "numpy · scipy · pandas", "click (required)"],
                  **MOD)],
            [Card(["matplotlib", "[plot] extra, optional", ""], **OPT),
             Card(["PySide6", "required by the Studio", "distribution"], **OPT),
             Card(["external executables", "IQ-TREE · MEGA-CC · FigTree",
                   "installed separately"], **OPT)],
        ], 4),
    ]

    # (from lane, from row, col, to lane, to row, col, kind); every pair is
    # column-aligned, so each arrow is a short straight line onto a box centre.
    #   data -> thin solid (data flow)      dep -> thin dashed (dependency)
    # An arrow always leaves the module that *uses* and points at the module
    # that is *used*, so it reproduces one real ``from .X import`` line (or the
    # value passed at run time, for the single data-flow arrow).  In particular
    # ``treeio.py:12 from .taxonomy import resolve_taxa_tokens`` runs treeio ->
    # taxonomy: taxonomy is never fed by treeio's callers, and ``table.py``
    # imports only ``report`` and ``tree``, so no arrow touches treeio from the
    # optional lane.
    EDGES = [
        (2, 0, 0, 1, 0, 0, "data"),         # treeio    -> rrf       (PhyloNode)
        (1, 0, 2, 2, 0, 2, "dep"),          # ci        -> report    (result types)
        (2, 0, 0, 3, 0, 0, "dep"),          # treeio    -> taxonomy  (resolve_taxa_tokens)
        (3, 0, 1, 2, 0, 1, "dep"),          # megacc    -> _external
        (3, 0, 2, 2, 0, 2, "dep"),          # viz       -> report
        (3, 1, 0, 4, 0, 0, "dep"),          # ddbd      -> tree
        (3, 1, 1, 4, 0, 1, "dep"),          # corrtest  -> _constants
        (3, 1, 2, 4, 0, 2, "dep"),          # bootstrap -> third-party stack
    ]
    KEY_LINES = ["thick arrow = RRF \u2192 calibration \u2192 interval pipeline   ·   "
                 "thin arrow = data flow",
                 "dashed arrow = dependency (uses)   ·   dashed box and grey fill = "
                 "optional module or extra"]

    def _wrap_lane_title(title):
        words = title.split()
        lines, cur = [], ""
        for wd in words:
            trial = (cur + " " + wd).strip()
            if text_mm(trial, M.FONT_HEAD, True) > LBL_W - 2 * CELL_PAD and cur:
                lines.append(cur)
                cur = wd
            else:
                cur = trial
        if cur:
            lines.append(cur)
        return lines

    def _geometry():
        """Lane bands, the dispatch-bus line and the total page height."""
        hs = [LANE_PAD * 2 + len(rows) * BOX_H + (len(rows) - 1) * ROW_GAP
              for _, rows, _ in LANES]
        ys, y = [], FRAME_INSET + 2.4
        for i, h in enumerate(hs):
            ys.append(y)
            y += h + LANE_GAP
            if i == 0:
                y += BUS_GAP - LANE_GAP          # room for the dispatch bus
        bus_y = ys[0] + hs[0] + (ys[1] - ys[0] - hs[0]) / 2.0
        key_y = y - LANE_GAP + 2.2
        page_h = key_y + 2.6 + len(KEY_LINES) * LINE + PAD + FRAME_INSET
        return hs, ys, bus_y, key_y, page_h

    hs, ys, bus_y, key_y, page_h = _geometry()
    sheet = Sheet(page_h, "Figure S1")

    # outer frame: >= 3 mm from every page edge
    sheet.box(FX0, FRAME_INSET, FX1 - FX0, page_h - 2 * FRAME_INSET,
              fill="white", edge=GREY_TXT, lw=LINE_MIN + 0.1)

    def top(l, r):
        return ys[l] + LANE_PAD + r * (BOX_H + ROW_GAP)

    def cx(c):
        return COL_X[c] + COL_W / 2.0

    for li, (title, rows, ramp) in enumerate(LANES):
        h = hs[li]
        sheet.box(FX0 + 1.6, ys[li], LBL_W, h, fill=RAMP[ramp],
                  edge=RAMP_EDGE[ramp], lw=LINE_MIN + 0.1)
        sheet.lines(FX0 + 1.6, ys[li], LBL_W, _wrap_lane_title(title),
                    pt=M.FONT_HEAD, bold=True, pad=CELL_PAD)
        if title.startswith("Optional"):         # dashed grouping rectangle
            gx = FX0 + 1.6 + LBL_W + 1.4
            sheet.box(gx, ys[li] + 1.0, AREA_X1 - gx - 0.2, h - 2.0,
                      fill="none", edge=GREY_TXT, lw=0.9, dashed=True)
        for r, row in enumerate(rows):
            for c, card in enumerate(row):
                lines = [l for l in card.lines if l]
                assert len(lines) <= N_LINES, lines
                card.lines = lines + [""] * (N_LINES - len(lines))
                card.height = BOX_H              # equal box heights
                card.draw(sheet, COL_X[c], top(li, r), COL_W)

    # ---- interfaces -> dispatch bus -> the three core engines (col centres) --
    sheet.line([(cx(0), bus_y), (cx(2), bus_y)], lw=W_DEP)
    for c in range(3):
        sheet.arrow([(cx(c), top(0, 0) + BOX_H), (cx(c), bus_y)], lw=W_DEP)
        sheet.arrow([(cx(c), bus_y), (cx(c), top(1, 0) - 0.35)], lw=W_DEP)
    # ---- core pipeline, centre to centre ---------------------------------- #
    for a in range(2):
        ymid = top(1, 0) + BOX_H / 2.0
        sheet.arrow([(COL_X[a] + COL_W + 0.4, ymid), (COL_X[a + 1] - 0.4, ymid)],
                    lw=W_STAGE)
    # ---- inter-lane semantic edges ----------------------------------------- #
    for l1, r1, c1, l2, r2, c2, kind in EDGES:
        ls = "-" if kind == "data" else DASH_STYLE
        if l2 > l1:
            a, b = (cx(c1), top(l1, r1) + BOX_H), (cx(c2), top(l2, r2) - 0.35)
        else:
            a, b = (cx(c1), top(l1, r1)), (cx(c2), top(l2, r2) + BOX_H + 0.35)
        assert abs(a[0] - b[0]) < 0.01, "edges stay on their column centre"
        sheet.arrow([a, b], lw=W_DEP, ls=ls)

    # ---- line-weight / shape grammar key ----------------------------------- #
    sheet.line([(FX0 + 1.6, key_y), (FX1 - 1.6, key_y)], lw=LINE_MIN,
               color=GREY_TXT)
    for i, kl in enumerate(KEY_LINES):
        sheet.text(FX0 + 1.6 + CELL_PAD, key_y + 2.6 + i * LINE, kl,
                   color=GREY_TXT, room=FX1 - FX0 - 3.2 - 2 * CELL_PAD)

    out_dir = HERE if out_dir is None else out_dir
    stem = os.path.join(out_dir, "附件图S1-软件架构图")
    pdf, png = sheet.save(stem, "OpenRelTime Figure S1 - software architecture",
                          "OpenRelTime; software architecture; modules; "
                          "Bioinformatics; supplementary figure")
    write_drawio(sheet, stem + ".drawio", "OpenRelTime-Architecture")
    rep = verify(pdf, sheet)
    rep["page_h_mm"] = round(page_h, 2)
    rep["frame_inset_mm"] = FRAME_INSET
    rep["lane_heights_mm"] = [round(h, 1) for h in hs]
    return rep, M


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--figure-engine", default=None,
                    help="folder holding make_figures.py (shared with Figure 1);"
                         " also read from $%s" % ENV_VAR)
    ap.add_argument("--outdir", default=None,
                    help="where to write the .drawio/.pdf/.png (default: here)")
    args = ap.parse_args()
    rep, M = build(args.outdir, args.figure_engine)
    if M.FAILS:
        print("%d fit or delivery failures:" % len(M.FAILS))
        for f in M.FAILS:
            print("  " + f)
    if os.environ.get("OPENRELTIME_FIT_COLLECT"):
        M.dump_measures()
    if M.FAILS:
        return 1
    print("Figure S1 (附件图S1-软件架构图) measured from the delivered PDF/PNG:")
    print(M._fmt(rep)
          + "\n  lane frame inset         {:.2f} mm  (limit >= 3.00)"
            "\n  lane heights (mm)        {}".format(
                rep["frame_inset_mm"], rep["lane_heights_mm"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
