# -*- coding: utf-8 -*-
"""OpenRelTime Figure 1 (analysis workflow) and the shared figure engine.

This module is the single source of truth for the two grammar-style figures of
the manuscript, Figure 1 (``正文图1-分析流程图``) and Figure S1 (drawn by
``../OpenRelTime-手稿附图/附件图S1-软件架构图.py``, which imports this engine):

  * the figure content is held as an explicit Python spec (``FIG1``, below) --
    every box, text run and arrow with its position in millimetres;
  * :class:`Sheet` renders that spec with matplotlib into a PDF (embedded
    TrueType fonts, so the file carries a real, searchable text layer) and a
    PNG at 1200 dpi, both exactly ``PAGE_W`` = 178 mm wide, the double-column
    width of the journal page;
  * :func:`write_drawio` mirrors the same spec into an editable ``.drawio``
    file, so an editor can move a box without touching the script;
  * :func:`verify` measures the *delivered* files (page box, raster width and
    ``pHYs`` density, ``/Font`` objects, extractable characters) and re-checks
    that no artwork crosses the page boundary, failing loudly if it does.

Running the script regenerates Figure 1 in place:

    python make_figures.py                  # drawio + PDF + PNG for Figure 1
    python make_figures.py --outdir /tmp    # write the three files elsewhere
    OPENRELTIME_FIT_COLLECT=1 python make_figures.py   # list every text run
        # with its measured width and the room available, for re-flowing copy.

Typography follows the manuscript: Unicode >= (U+2265), congruent-to, en dashes
for ranges and parenthetical dashes, comma thousands separators and ``10%``
(digit + percent sign, no space).  Capability comparison lives in Supplementary
Table S3 only; it is deliberately not repeated inside this figure.
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, PathPatch  # noqa: E402
from matplotlib.path import Path  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# --------------------------------------------------------------------------- #
# printing constraints
# --------------------------------------------------------------------------- #
MM = 1.0 / 25.4
PAGE_W = 178.0                 # mm, double-column width (504.57 pt)
PX_PER_MM = 96.0 / 25.4        # drawio page units: 178 mm -> 673 px
FRAME_INSET = 3.6              # mm, minimum distance of any frame from a page edge
LINE = 3.2275                  # mm, baseline pitch of a FONT_BODY line: the pitch
                               # of the approved artwork, which also puts every
                               # card height on the 96 dpi grid of the drawio mirror
PAD = 1.5                      # mm, inner padding of a card (vertical)
INSET = 1.6                    # mm, inner padding of a card (horizontal)
PNG_DPI = 1200.0
MIN_PT = 7.5                   # smallest type size allowed in a delivered figure

FONT_BODY = 7.5                # pt, card text
FONT_HEAD = 8.5                # pt, lane / stage titles

# line widths in points; the drawio mirror converts pt -> px at 96/72
LINE_MIN = 0.75
LINE_CARD = 0.85
LINE_OPT = 0.90
LINE_HEAD = 0.95
LINE_THICK = 1.50
W_DEP = 1.00                   # thin arrow: data flow / dependency
W_STAGE = 1.50                 # thick arrow: stage order

DASH_STYLE = (0, (4.0, 2.8))   # optional module / opt-in path
DOT_STYLE = (0, (0.1, 1.6))    # external tool call

INK = "#111111"
GREY_TXT = "#333333"
CARD_EDGE = "#111111"
DOC_FILL = "#F5F2E8"
DOC_EDGE = "#5B7C99"
OPT_FILL = "#F2F2F2"
LEG_FILL = "#F7F7F7"
RAMP = ["#DDE7EF", "#CCD9E5", "#BACADF", "#A8BBDA", "#96ACD0"]
RAMP_EDGE = ["#5B7C99", "#4E6F8C", "#41627F", "#345572", "#274865"]

matplotlib.rcParams.update({
    "pdf.fonttype": 42,        # embedded TrueType -> searchable text layer
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.size": FONT_BODY,
    "savefig.dpi": PNG_DPI,
    "figure.dpi": 100,
})

# --------------------------------------------------------------------------- #
# text metrics and the fit collector
# --------------------------------------------------------------------------- #
COLLECT = [bool(os.environ.get("OPENRELTIME_FIT_COLLECT"))]
FAILS: list = []
_MEASURES: list = []

_probe_fig = None
_probe_txt = None
_probe_rend = None


def _probe():
    global _probe_fig, _probe_txt, _probe_rend
    if _probe_fig is None:
        _probe_fig = plt.figure(figsize=(6, 6), dpi=72.0)
        _probe_fig.add_axes([0, 0, 1, 1]).axis("off")
        _probe_rend = _probe_fig.canvas.get_renderer()
        _probe_txt = _probe_fig.text(0, 0, "")
    return _probe_txt, _probe_rend


def text_mm(s, pt=FONT_BODY, bold=False):
    """Width of a single text run in millimetres, as matplotlib will draw it."""
    t, r = _probe()
    t.set_text(s)
    t.set_fontsize(pt)
    t.set_fontweight("bold" if bold else "normal")
    return float(t.get_window_extent(r).width) / 72.0 * 25.4


def line_height_mm(pt=FONT_BODY):
    """Height of one text line in millimetres (matplotlib default spacing)."""
    t, r = _probe()
    t.set_text("Hp")
    t.set_fontsize(pt)
    return float(t.get_window_extent(r).height) / 72.0 * 25.4


def _fit(run, w, room, where):
    """Record a text run; complain when it does not fit the room it was given."""
    _MEASURES.append((where, run, round(w, 2), None if room is None else round(room, 2)))
    if room is not None and w > room + 0.02:
        FAILS.append("text overflows its room by %.2f mm at %s: %r"
                     % (w - room, where, run[:56]))


# --------------------------------------------------------------------------- #
# cards
# --------------------------------------------------------------------------- #
class Card:
    """A labelled box: one border plus its text lines.

    ``lines`` entries are plain strings (first line bold when ``bold_first``) or
    ``(text, bold)`` tuples.  Empty strings are spacers and draw nothing, which
    lets a caller keep every card on one grid of equal height.
    """

    def __init__(self, lines, height=0.0, group="main", bold_first=True,
                 shape="box", fill="white", edge=None, lw=None, dashed=None,
                 pt=FONT_BODY, pitch=None, inset=INSET, pad=None, id=None):
        self.lines = list(lines)
        self.height = height
        self.group = group
        self.bold_first = bold_first
        self.shape = shape
        self.fill = fill
        self.edge = edge
        self.lw = lw
        self.dashed = dashed
        self.pt = pt
        self.pitch = LINE if pitch is None else pitch
        self.inset = inset
        self.pad = PAD if pad is None else pad
        self.id = id

    def draw(self, sheet, x, y, w):
        dashed = (self.group == "optional") if self.dashed is None else self.dashed
        lw = self.lw
        if lw is None:
            lw = LINE_OPT if dashed else LINE_CARD
        edge = self.edge
        if edge is None:
            edge = DOC_EDGE if self.shape == "document" else CARD_EDGE
        sheet.box(x, y, w, self.height, fill=self.fill, edge=edge, lw=lw,
                  dashed=dashed, shape=self.shape, id=self.id)
        for i, ln in enumerate(self.lines):
            run, bold = (ln, self.bold_first and i == 0) if isinstance(ln, str) else ln
            # an empty line is still recorded: it keeps the card's text block
            # aligned with the grid and shows up in the drawio mirror
            sheet.text(x + self.inset, y + self.pad + i * self.pitch, run,
                       pt=self.pt, bold=bold, room=(None if not run
                                                    else w - 2 * self.inset))


# --------------------------------------------------------------------------- #
# the sheet
# --------------------------------------------------------------------------- #
class Sheet:
    """A page of ``PAGE_W`` x ``page_h`` millimetres, origin at the top left."""

    def __init__(self, page_h, name="Figure"):
        self.page_h = float(page_h)
        self.name = name
        self.items = []
        self.min_pt = None
        self.fig = plt.figure(figsize=(PAGE_W * MM, self.page_h * MM), dpi=100)
        self.ax = self.fig.add_axes([0.0, 0.0, 1.0, 1.0])
        self.ax.set_xlim(0.0, PAGE_W)
        self.ax.set_ylim(self.page_h, 0.0)
        self.ax.axis("off")
        self.fig.patch.set_facecolor("white")

    # -- primitives -------------------------------------------------------- #
    def box(self, x, y, w, h, fill="white", edge=INK, lw=LINE_CARD, dashed=False,
            shape="box", id=None):
        it = dict(kind="box", x=x, y=y, w=w, h=h, fill=fill, edge=edge, lw=lw,
                  dashed=dashed, shape=shape)
        if id is not None:
            it["id"] = id
        self.items.append(it)
        r = 0.55
        if shape == "document":
            e = min(2.9, 0.45 * w, 0.6 * h)          # the dog-ear
            body = [(x, y), (x + w - e, y), (x + w, y + e), (x + w, y + h),
                    (x, y + h), (x, y)]
            it["artist"] = self.ax.add_patch(PathPatch(
                Path(body, [Path.MOVETO] + [Path.LINETO] * 4 + [Path.CLOSEPOLY]),
                facecolor=fill, edgecolor=edge, lw=lw,
                linestyle=DASH_STYLE if dashed else "-", zorder=2))
            self.ax.add_patch(PathPatch(
                Path([(x + w - e, y), (x + w - e, y + e), (x + w, y + e),
                      (x + w - e, y + e)],
                     [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO]),
                facecolor="none", edgecolor=edge, lw=lw,
                linestyle=DASH_STYLE if dashed else "-", zorder=3))
        else:
            it["artist"] = self.ax.add_patch(FancyBboxPatch(
                (x, y), w, h, boxstyle="round,pad=0,rounding_size=%g" % r,
                facecolor=fill, edgecolor=edge, lw=lw,
                linestyle=DASH_STYLE if dashed else "-",
                mutation_aspect=1.0, zorder=2))
        return it

    def text(self, x, y, s, pt=FONT_BODY, bold=False, color=INK, room=None,
             id=None):
        w = text_mm(s, pt, bold)
        h = line_height_mm(pt)
        _fit(s, w, room, "y=%.2f" % y)
        it = dict(kind="text", x=x, y=y, s=s, pt=pt, bold=bold, color=color,
                  w=w, h=h)
        if id is not None:
            it["id"] = id
        self.items.append(it)
        it["artist"] = self.ax.text(x, y, s, fontsize=pt, color=color,
                                    fontweight="bold" if bold else "normal",
                                    ha="left", va="top", zorder=5)
        self.min_pt = pt if self.min_pt is None else min(self.min_pt, pt)
        return it

    def lines(self, x, y, w, ls, pt=FONT_BODY, bold=False, pad=INSET, room=None):
        """A block of lines inside a box whose left edge and width are given."""
        avail = w - 2 * pad if room is None else room
        out = []
        for i, ln in enumerate(ls):
            run, b = (ln, bold) if isinstance(ln, str) else (ln[0], ln[1] or bold)
            if not run:
                continue
            out.append(self.text(x + pad, y + pad + i * (pt * 0.3528 * 1.2), run,
                                 pt=pt, bold=b, room=avail))
        return out

    def line(self, pts, lw=W_DEP, color=INK, ls="-", id=None):
        it = dict(kind="line", pts=[tuple(p) for p in pts], lw=lw, color=color,
                  ls=ls, arrow=False)
        if id is not None:
            it["id"] = id
        self.items.append(it)
        it["artist"] = self._stroke(pts, lw, color, ls)
        return it

    def arrow(self, pts, lw=W_DEP, color=INK, ls="-", id=None):
        pts = [tuple(p) for p in pts]
        it = dict(kind="arrow", pts=pts, lw=lw, color=color, ls=ls, arrow=True)
        if id is not None:
            it["id"] = id
        self.items.append(it)
        it["artist"] = self._stroke(pts, lw, color, ls)
        (x0, y0), (x1, y1) = pts[-2], pts[-1]
        self.ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=4,
                         arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                         linestyle=ls, shrinkA=0, shrinkB=0,
                                         mutation_scale=9.0))
        return it

    def _stroke(self, pts, lw, color, ls):
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        return self.ax.plot(xs, ys, lw=lw, color=color, ls=ls, solid_capstyle="butt",
                            zorder=3)[0]

    # -- output ------------------------------------------------------------ #
    def save(self, stem, title, keywords):
        pdf, png = stem + ".pdf", stem + ".png"
        self.fig.canvas.draw()
        self.fig.savefig(png, dpi=PNG_DPI, facecolor="white")
        self.fig.savefig(pdf, facecolor="white")
        set_pdf_metadata(pdf, title, keywords)
        return pdf, png

    # -- self-checks ------------------------------------------------------- #
    def clipped(self, tol_px=0.6):
        """Artists whose bounding box crosses the page rectangle."""
        self.fig.canvas.draw()
        rend = self.fig.canvas.get_renderer()
        w_px, h_px = self.fig.canvas.get_width_height()
        bad = []
        for it in self.items:
            a = it.get("artist")
            if a is None:
                continue
            try:
                bb = a.get_window_extent(rend)
            except TypeError:
                bb = a.get_window_extent()
            if bb.x0 < -tol_px or bb.y0 < -tol_px or bb.x1 > w_px + tol_px \
                    or bb.y1 > h_px + tol_px:
                bad.append((it["kind"], repr(_label(it))[:48],
                            round(bb.x0, 1), round(bb.y0, 1),
                            round(bb.x1, 1), round(bb.y1, 1)))
        return bad, (w_px, h_px)


def _label(it):
    return it.get("s") or it.get("shape") or it["kind"]


# --------------------------------------------------------------------------- #
# drawio mirror
# --------------------------------------------------------------------------- #
def _px(v):
    return "%.1f" % round(v * PX_PER_MM, 1)


def _num(v):
    s = "%.2f" % v
    return s[:-1] if s.endswith("0") and "." in s else s


def _esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
             .replace('"', "&quot;"))


def _inside(box, txt):
    return (box["x"] - 0.05 <= txt["x"]
            and txt["x"] + txt["w"] <= box["x"] + box["w"] + 0.05
            and box["y"] - 0.05 <= txt["y"]
            and txt["y"] + txt["h"] <= box["y"] + box["h"] + 0.05)


BOX_STYLE = ("rounded=1;html=1;whiteSpace=wrap;fillColor=%(fill)s;strokeColor=%(edge)s;"
             "strokeWidth=%(lw)s;arcSize=10;fontColor=#111111;align=left;"
             "verticalAlign=top;spacingLeft=6;spacingTop=6;spacingBottom=9;fontSize=10;")
DOC_STYLE = BOX_STYLE.replace("rounded=1;", "shape=document;")
OPT_STYLE = "dashed=1;dashPattern=6 4;"
FREE_STYLE = ("text;html=1;align=left;verticalAlign=middle;fontColor=%(color)s;"
              "fontSize=10;")
EDGE_STYLE = ("endArrow=%(head)s;html=1;strokeColor=%(color)s;strokeWidth=%(lw)s;"
              "rounded=0;")


def write_drawio(sheet, path, name):
    """Emit an editable draw.io mirror of the sheet (same page, same geometry).

    Text runs are folded into every box that contains them, sorted top-to-bottom
    and left-to-right, which is how the figure reads in the draw.io editor; runs
    that belong to no box become free text cells.
    """
    boxes = [it for it in sheet.items if it["kind"] == "box"]
    texts = [it for it in sheet.items if it["kind"] == "text"]
    edges = [it for it in sheet.items if it["kind"] in ("line", "arrow")]
    for b in boxes:
        # runs keep the order in which they were drawn, which is how the card
        # they belong to wrote them: line by line, block by block
        b["runs"] = [t for t in texts if _inside(b, t)]
    free = [t for t in texts if not any(_inside(b, t) for b in boxes)]

    used = [it["id"] for it in boxes + free + edges if "id" in it]
    nxt = (max(used) + 1) if used else 101
    for it in boxes + free + edges:
        if "id" not in it:
            it["id"] = nxt
            nxt += 1

    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<!-- GENERATED by make_figures.py from the same spec that renders the',
           '     PDF and the PNG.  Edit the script and re-run it; do not hand-edit',
           '     this file.  Page = 178 mm at 96 dpi, so printed sizes match the PDF. -->',
           '<mxfile host="drawio" version="26.0.0">',
           '  <diagram name="%s">' % name,
           '    <mxGraphModel dx="0" dy="0" grid="1" gridSize="10" guides="1"'
           ' tooltips="1" connect="1" arrows="1" fold="1" page="1" pageScale="1"'
           ' pageWidth="%d" pageHeight="%d" math="0" shadow="0" background="#ffffff">'
           % (int(round(PAGE_W * PX_PER_MM)), int(round(sheet.page_h * PX_PER_MM))),
           '      <root>',
           '        <mxCell id="0" />',
           '        <mxCell id="1" parent="0" />']

    def cell(it, style, value):
        g = it["geometry"]
        out.append('        <mxCell id="%s" value="%s" style="%s" vertex="1"'
                   ' parent="1">' % (it["id"], _esc(value), style))
        out.append('          %s' % g)
        out.append('        </mxCell>')

    for b in boxes:
        st = (DOC_STYLE if b["shape"] == "document" else BOX_STYLE) % {
            "fill": b["fill"], "edge": b["edge"], "lw": _num(b["lw"] * 96.0 / 72.0)}
        if b["dashed"]:
            st += OPT_STYLE
        b["geometry"] = ('<mxGeometry x="%s" y="%s" width="%s" height="%s"'
                         ' as="geometry" />' % (_px(b["x"]), _px(b["y"]),
                                                _px(b["w"]), _px(b["h"])))
        cell(b, st, "\n".join(t["s"] for t in b["runs"]))
    for t in free:
        t["geometry"] = ('<mxGeometry x="%s" y="%s" width="%s" height="%s"'
                         ' as="geometry" />' % (_px(t["x"]), _px(t["y"]),
                                                _px(t["w"]), _px(t["h"])))
        cell(t, FREE_STYLE % {"color": t["color"]}, t["s"])
    for e in edges:
        pts = e["pts"]
        st = EDGE_STYLE % {"head": "classic" if e["arrow"] else "none",
                           "color": e["color"], "lw": _num(e["lw"] * 96.0 / 72.0)}
        if e["ls"] != "-":
            st += OPT_STYLE
        seg = ['<mxPoint x="%s" y="%s" as="sourcePoint" />' % (_px(pts[0][0]), _px(pts[0][1]))]
        for p in pts[1:-1]:
            seg.append('<mxPoint x="%s" y="%s" />' % (_px(p[0]), _px(p[1])))
        seg.append('<mxPoint x="%s" y="%s" as="targetPoint" />'
                   % (_px(pts[-1][0]), _px(pts[-1][1])))
        out.append('        <mxCell id="%s" style="%s" edge="1" parent="1">'
                   % (e["id"], st))
        out.append('          <mxGeometry relative="1" as="geometry">%s</mxGeometry>'
                   % "".join(seg))
        out.append('        </mxCell>')
    out += ['      </root>', '    </mxGraphModel>', '  </diagram>', '</mxfile>', '']
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))
    return path


# --------------------------------------------------------------------------- #
# delivery checks
# --------------------------------------------------------------------------- #
def set_pdf_metadata(pdf_path, title, keywords):
    """Readable DocInfo/XMP metadata on the exported PDF (needs PyMuPDF)."""
    try:
        import fitz
    except ImportError:
        return False
    doc = fitz.open(pdf_path)
    md = doc.metadata or {}
    md.update({"title": title, "author": "", "creator": "OpenRelTime",
               "producer": "matplotlib", "keywords": keywords})
    doc.set_metadata(md)
    tmp = pdf_path + ".tmp"
    doc.save(tmp)
    doc.close()
    os.replace(tmp, pdf_path)
    return True


def verify(pdf, sheet):
    """Measure the delivered PDF and PNG; register a FAIL for each breach."""
    import fitz
    from PIL import Image

    rep = {}
    doc = fitz.open(pdf)
    page = doc[0]
    rep["pdf_w_mm"] = page.rect.width * 25.4 / 72.0
    rep["pdf_h_mm"] = page.rect.height * 25.4 / 72.0
    rep["pdf_pt"] = "%.2f x %.2f" % (page.rect.width, page.rect.height)
    with open(pdf, "rb") as fh:
        rep["fonts"] = fh.read().count(b"/Font")
    rep["chars"] = len(page.get_text().strip())
    doc.close()

    png = pdf[:-4] + ".png"
    im = Image.open(png)
    rep["png_w"], rep["png_h"] = im.size
    info = im.info
    dpi = info.get("dpi") or info.get("physical dpi_x")
    rep["png_dpi"] = float(dpi[0]) if isinstance(dpi, (tuple, list)) else float(dpi or 0.0)
    im.close()

    bad, canvas = sheet.clipped()
    rep["clipped"] = bad
    rep["min_pt"] = sheet.min_pt

    if abs(rep["pdf_w_mm"] - PAGE_W) > 0.02:
        FAILS.append("PDF width %.3f mm != %.2f mm" % (rep["pdf_w_mm"], PAGE_W))
    if abs(rep["pdf_h_mm"] - sheet.page_h) > 0.02:
        FAILS.append("PDF height %.3f mm != spec %.3f mm" % (rep["pdf_h_mm"], sheet.page_h))
    if rep["fonts"] < 1:
        FAILS.append("PDF has no /Font object: the text layer is missing")
    if rep["chars"] < 40:
        FAILS.append("PDF text layer extracts only %d characters" % rep["chars"])
    want_px = int(round(PAGE_W / 25.4 * PNG_DPI))
    if rep["png_w"] != want_px:
        FAILS.append("PNG is %d px wide, expected %d" % (rep["png_w"], want_px))
    if rep["png_dpi"] < PNG_DPI - 1.0:
        FAILS.append("PNG pHYs declares %.1f dpi, expected %.0f"
                     % (rep["png_dpi"], PNG_DPI))
    if bad:
        FAILS.append("%d element(s) cross the page boundary" % len(bad))
    if rep["min_pt"] is not None and rep["min_pt"] < MIN_PT - 1e-9:
        FAILS.append("smallest type %.2f pt is below the %.2f pt limit"
                     % (rep["min_pt"], MIN_PT))
    return rep


def _fmt(rep):
    return ("  page size (PDF)          %s pt = %.2f x %.2f mm  (width limit %.2f)\n"
            "  PNG raster               %d x %d px  at %.1f dpi (pHYs present)\n"
            "  PDF text layer           %d characters, %d /Font objects\n"
            "  page-bound check         %d element(s) clipped at the page edge\n"
            "  smallest type            %.2f pt (limit %.2f)"
            % (rep["pdf_pt"], rep["pdf_w_mm"], rep["pdf_h_mm"], PAGE_W,
               rep["png_w"], rep["png_h"], rep["png_dpi"],
               rep["chars"], rep["fonts"], len(rep["clipped"]),
               rep["min_pt"] or 0.0, MIN_PT))


def dump_measures():
    print("\n  text runs (measured width vs room):")
    for where, run, w, room in sorted(_MEASURES, key=lambda m: -(m[2] - (m[3] or 1e9))):
        print("    %6.2f / %s  %s" % (w, ("%.2f" % room) if room else "  -  ",
                                      run[:74]))


# --------------------------------------------------------------------------- #
# Figure 1 -- the spec
# --------------------------------------------------------------------------- #
# Geometry is in millimetres from the top-left page corner.  Column grid:
COL = [3.2, 61.59, 120.0]
CW = 54.8                       # card width; a card's text room is CW - 2*INSET
STAGE_H = 8.0
SHORT_H = 9.45              # a two-line card
TALL_H = 16.27
MID_H = 12.67
OPT_H = 13.04
NOTE_X = 63.1                   # the side note between the two input files
NOTE_DY = 3.23                  # its line pitch

# The bottom band: the capability table that used to sit here is gone on purpose
# (it duplicated Supplementary Table S3, contradicted it in one cell and pushed
# the figure to 185 mm); the external-tool box now spans the full width.
# The band carries two groups, because the engine treats the three programs
# differently (_external.py:16-19 registers only IQ-TREE for `bootstrap` and
# MEGA-CC for `megacc`, and `subprocess` appears nowhere else in the package):
# IQ-TREE and MEGA-CC are started as subprocesses from a validated argument
# list, while FigTree is only the downstream viewer whose NEXUS format the
# writers target and is never invoked.  Hence two headings on one title row,
# split by a rule, instead of a single "run as subprocesses" label.
BAND_Y = 141.0
BAND_H = 16.7
BAND_TITLE_Y = 141.75             # umbrella heading
BAND_GROUP_Y = 144.85             # the two role headings, one per group
BAND_CARD_Y = 148.0
BAND_CARD_H = 8.6
BAND_CARDS = [4.8, 61.2, 117.6]
BAND_SPLIT_X = 116.8              # rule between the two groups
PAGE_H_1 = BAND_Y + BAND_H + 3.2


def _stage(sheet, cid, ramp, x, y, w, num, label):
    """A stage header: the number and its name on one centred, bold line."""
    sheet.box(x, y, w, STAGE_H, fill=RAMP[ramp], edge=RAMP_EDGE[ramp], lw=LINE_HEAD,
              id=cid)
    ty = y + (STAGE_H - line_height_mm(FONT_HEAD)) / 2.0
    sheet.text(x + INSET, ty, num, pt=FONT_HEAD, bold=True, room=6.0)
    sheet.text(x + INSET + text_mm(num + " ", FONT_HEAD, True), ty, label,
               pt=FONT_HEAD, bold=True, room=w - 2 * INSET - text_mm(num + " ", FONT_HEAD, True))


def _fig1(sheet):
    """Draw Figure 1 from its spec, in the order the drawio mirror numbers cells."""
    # -- stage 1: input and validation ------------------------------------- #
    # The two input cards keep the wording of the audited artwork, re-flowed so
    # that every line fits its card in DejaVu Sans (the font Figures S1-S3 use);
    # the previous render relied on draw.io's narrower Helvetica metrics.
    Card(["Branch-length tree, Newick or",
          "NEXUS; substitutions per site"], SHORT_H, shape="document", bold_first=False,
         fill=DOC_FILL, edge=DOC_EDGE, lw=LINE_CARD, id=101).draw(sheet, COL[0], 3.2, CW)
    Card(["Calibration file (TSV): min/max",
          "bounds or probability densities"], SHORT_H, shape="document", bold_first=False,
         fill=DOC_FILL, edge=DOC_EDGE, lw=LINE_CARD, id=102).draw(sheet, COL[2], 3.2, CW)
    _stage(sheet, 103, 0, COL[0], 16.67, CW, "1", "Input & validation")
    Card(["CORE",
          "Root on the outgroup MRCA,",
          "prune it; require a binary",
          "ingroup, lengths \u2265 0"], TALL_H, lw=LINE_THICK, id=104).draw(sheet, COL[0], 27.07, CW)
    Card(["OPTIONAL",
          "Polytomies: refused by default",
          "(opt-in random binarisation",
          "on zero-length stems)"], TALL_H, group="optional", id=105).draw(sheet, COL[0], 45.03, CW)
    Card(["Zero-length branches: warn at",
          "\u2265 10% of branches (R3F)"], SHORT_H, group="optional",
         bold_first=False, id=106).draw(sheet, COL[0], 63.02, CW)

    # -- stage 2: the RRF core --------------------------------------------- #
    _stage(sheet, 107, 1, COL[1], 16.67, CW, "2", "RRF core")
    Card(["CORE",
          "Postorder: closed-form local",
          "rate at every node; subtrees",
          "folded into composite lengths"], TALL_H, lw=LINE_THICK, id=108).draw(sheet, COL[1], 27.07, CW)
    Card(["Preorder: multiply by the",
          "ancestral rate, divide times",
          "(grandparent rule)"], MID_H, lw=LINE_THICK, bold_first=False,
         id=109).draw(sheet, COL[1], 45.03, CW)
    Card(["Relative rates and node times",
          "(ingroup root rate \u2261 1)"], SHORT_H, lw=LINE_THICK, bold_first=False,
         id=110).draw(sheet, COL[1], 59.43, CW)

    # -- stages 3 and 4: calibration, intervals ---------------------------- #
    _stage(sheet, 111, 2, COL[2], 16.67, CW, "3", "Calibration")
    Card(["CORE",
          "Global time factor f = midpoint",
          "of the feasible [min, max]",
          "intersection"], TALL_H, lw=LINE_THICK, id=112).draw(sheet, COL[2], 27.07, CW)
    Card(["Violations: rescale and",
          "propagate the lineage rate,",
          "iterate (cap 100) or report"], MID_H, lw=LINE_THICK, bold_first=False,
         id=113).draw(sheet, COL[2], 45.03, CW)
    Card(["OPTIONAL",
          "Density calibrations: effective",
          "bounds, 10,000 replicates"], OPT_H, group="optional", id=114).draw(sheet, COL[2], 59.43, CW)
    _stage(sheet, 115, 3, COL[2], 75.27, CW, "4", "Intervals & options")
    Card(["CORE",
          "Analytical confidence intervals",
          "(delta method), truncated at",
          "the hard calibration bounds"], TALL_H, lw=LINE_THICK, id=116).draw(sheet, COL[2], 85.67, CW)
    # 117 sits on the column edge (453.5 px), not 2.3 px inside it
    Card(["Independent option modules:",
          "no order, no interdependency",
          "\u2022 CorrTest: rate autocorrelation",
          "\u2022 ddBD: birth-death tree prior",
          "\u2022 Bag-of-little-bootstraps dates",
          "\u2022 MEGA-CC bridge (MEGA users)",
          "\u2022 viz: timetree, interval plots",
          "\u2022 table: tree2table records",
          "\u2022 taxonomy: monophyly check"], 32.86, group="optional",
         edge=GREY_TXT, id=117).draw(sheet, COL[2], 103.66, CW)

    # -- stage 5: output ---------------------------------------------------- #
    _stage(sheet, 118, 4, COL[0], 75.27, 113.19, "5", "Output")
    Card(["Rates, times and",
          "calibrated ages:",
          "CSV tables"], MID_H, shape="document", bold_first=False, fill=DOC_FILL,
         edge=DOC_EDGE, lw=LINE_CARD, id=119).draw(sheet, 4.71, 85.67, 34.53)
    Card(["Timetree with",
          "[&rate=...]",
          "annotations"], MID_H, shape="document", bold_first=False, fill=DOC_FILL,
         edge=DOC_EDGE, lw=LINE_CARD, id=120).draw(sheet, 40.43, 85.67, 34.53)
    Card(["JSON report:",
          "parameters, seed,",
          "warnings"], MID_H, shape="document", bold_first=False, fill=DOC_FILL,
         edge=DOC_EDGE, lw=LINE_CARD, id=121).draw(sheet, 76.17, 85.67, 34.53)
    Card(["Three parallel products of one run;",
          "API and command line write identical files"], SHORT_H, group="optional",
         bold_first=False, id=122).draw(sheet, 4.71, 100.07, 110.21)

    # -- grammar key -------------------------------------------------------- #
    sheet.box(3.2, 112.32, 113.19, 19.92, fill=LEG_FILL, edge=GREY_TXT,
              lw=LINE_CARD, id=123)
    sheet.box(4.71, 113.83, 7.01, 3.23, fill="white", edge=CARD_EDGE, lw=LINE_THICK, id=124)
    sheet.box(4.71, 118.07, 7.01, 3.23, fill="white", edge=CARD_EDGE, lw=LINE_CARD, id=125)
    sheet.box(4.71, 122.27, 7.01, 3.23, fill="white", edge=CARD_EDGE, lw=LINE_OPT,
              dashed=True, id=126)
    sheet.box(4.71, 126.50, 7.01, 3.23, fill="white", edge=CARD_EDGE, lw=LINE_CARD,
              shape="document", id=127)
    key_rows = [("core engine (cf. Fig. S1)", "stage order, 1 to 5"),
                ("processing step", "real data dependency"),
                ("option module / opt-in behaviour", "side input / optional path"),
                ("data object (file), not a step", "external tool call")]
    for i, (left, right) in enumerate(key_rows):
        y = 113.83 + i * 4.24 + 0.3
        sheet.text(13.32, y, left, room=53.0)
        sheet.text(69.30, y, right, room=45.4)

    # -- external tools, full-width band ------------------------------------ #
    sheet.box(3.2, BAND_Y, 171.6, BAND_H, fill="white", edge=GREY_TXT, lw=LINE_OPT,
              dashed=True, id=129)
    sheet.text(4.8, BAND_TITLE_Y,
               "External programs \u2013 installed separately and never imported"
               " by the engine", bold=True, room=171.6 - 2 * INSET, id=None)
    # two role headings, one per group: the engine starts IQ-TREE and MEGA-CC as
    # subprocesses, and writes files for FigTree without ever starting it.
    sheet.text(4.8, BAND_GROUP_Y,
               "Called as subprocesses from a validated argument list (no shell)",
               bold=True, room=BAND_SPLIT_X - 4.8 - 1.2, id=None)
    sheet.text(BAND_CARDS[2], BAND_GROUP_Y, "Downstream viewer: never called",
               bold=True, room=3.2 + 171.6 - INSET - BAND_CARDS[2], id=None)
    sheet.line([(BAND_SPLIT_X, BAND_GROUP_Y - 0.6),
                (BAND_SPLIT_X, BAND_CARD_Y + BAND_CARD_H + 0.2)],
               lw=LINE_MIN, color=GREY_TXT)
    for cid, (name, what) in zip((130, 131, 132),
                                 [("IQ-TREE", "runs the bag-of-little-bootstraps"),
                                  ("MEGA-CC", "runs the MEGA-CC bridge"),
                                  ("FigTree", "opens the annotated NEXUS output")]):
        Card([name, what], BAND_CARD_H, group="optional", fill=OPT_FILL,
             edge=CARD_EDGE, lw=LINE_OPT, bold_first=True, id=cid).draw(
            sheet, BAND_CARDS[cid - 130], BAND_CARD_Y, CW)

    # -- side note between the two input files ------------------------------ #
    for nid, dy, run in ((133, 0.0, "Two files are read per run: the tree"),
                         (134, NOTE_DY, "drives stages 1\u20132; the calibration file"),
                         (135, 2 * NOTE_DY, "is a side input to stages 3 and 4.")):
        sheet.text(NOTE_X, 6.32 + dy, run, color=GREY_TXT, id=nid)

    # -- arrows ------------------------------------------------------------- #
    sheet.arrow([(22.38, 12.65), (22.38, 16.25)], lw=W_DEP, ls=DASH_STYLE, id=136)
    sheet.arrow([(139.18, 12.65), (139.18, 16.25)], lw=W_DEP, ls=DASH_STYLE, id=137)
    sheet.arrow([(58.39, 20.67), (61.20, 20.67)], lw=W_STAGE, id=138)
    sheet.arrow([(116.79, 20.67), (119.60, 20.67)], lw=W_STAGE, id=139)
    sheet.arrow([(139.18, 72.87), (139.18, 74.88)], lw=W_STAGE, id=140)
    sheet.arrow([(159.46, 72.79), (159.46, 74.98)], lw=W_DEP, ls=DASH_STYLE, id=141)
    sheet.arrow([(61.30, 115.44), (67.71, 115.44)], lw=W_STAGE, id=142)
    sheet.arrow([(61.30, 119.68), (67.71, 119.68)], lw=W_DEP, id=143)
    sheet.arrow([(61.30, 123.89), (67.71, 123.89)], lw=W_DEP, ls=DASH_STYLE, id=144)
    sheet.arrow([(61.30, 128.11), (67.71, 128.11)], lw=LINE_OPT, ls=DOT_STYLE, id=145)
    sheet.arrow([(119.60, 79.27), (116.79, 79.27)], lw=W_STAGE, id=146)
    # 148: the external-tool call.  It leaves 1.6 mm (6 px) below the dashed
    # bottom border of 117 instead of riding on it, then elbows across the
    # corridor under the key so that it lands on the two programs the engine
    # actually starts (the gap between the IQ-TREE and MEGA-CC cards), not on
    # FigTree, which the engine only writes files for.
    sheet.arrow([(147.40, 138.12), (147.40, 139.40), (60.40, 139.40),
                 (60.40, 140.81)], lw=LINE_OPT, ls=DASH_STYLE, id=148)


def build_figure_1(out_dir):
    sheet = Sheet(PAGE_H_1, "Figure 1")
    _fig1(sheet)
    stem = os.path.join(out_dir, "正文图1-分析流程图")
    pdf, png = sheet.save(
        stem, "OpenRelTime Figure 1 - analysis workflow",
        "OpenRelTime; relative rate framework; molecular dating; workflow; "
        "Bioinformatics; Application Note; figure")
    write_drawio(sheet, stem + ".drawio", "OpenRelTime-Workflow")
    rep = verify(pdf, sheet)
    rep["page_h_mm"] = round(PAGE_H_1, 2)
    rep["drawio"] = os.path.basename(stem) + ".drawio"
    return rep


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default=HERE,
                    help="where to write the .drawio/.pdf/.png (default: this folder)")
    args = ap.parse_args()
    rep = build_figure_1(args.outdir)
    if FAILS:
        print("%d fit or delivery failures:" % len(FAILS))
        for f in FAILS:
            print("  " + f)
        return 1
    print("Figure 1 (正文图1-分析流程图) measured from the delivered PDF/PNG:")
    print(_fmt(rep))
    print("  figure height            %.2f mm  (was 185.36 before the capability"
          " table was removed)" % rep["page_h_mm"])
    for kind, lab, x0, y0, x1, y1 in rep["clipped"]:
        print("  clipped:", kind, lab, x0, y0, x1, y1)
    if COLLECT[0]:
        dump_measures()
    return 0


if __name__ == "__main__":
    sys.exit(main())
