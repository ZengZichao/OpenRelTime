# Figure sources

Source code and data for the four figures of the OpenRelTime manuscript:
**Figure 1** (analysis workflow), **Figure S1** (software architecture),
**Figure S2** (agreement with R3F and wall-clock time) and **Figure S3**
(effective bounds, calibration sensitivity, and the inconsistency in the
shipped example calibration file).  The file names are the manuscript's own,
unchanged, so they stay traceable to the submitted artwork.

| File | Produces |
|---|---|
| `make_figures.py` | 正文图1-分析流程图 (Figure 1) — and the shared drawing engine that Figure S1 imports |
| `附件图S1-软件架构图.py` | 附件图S1-软件架构图 (Figure S1) |
| `附件图S2-性能对比图.py` + `附件图S2-性能对比图-data.json` | 附件图S2-性能对比图 (Figure S2) |
| `附件图S3-置信区间图.py` + `附件图S3-置信区间图-data.json` | 附件图S3-置信区间图 (Figure S3) |

The rendered `.pdf`/`.png`/`.drawio` artwork is submitted with the manuscript
and is deliberately **not** checked in here; each script re-measures its own
delivered files (page width, raster dpi, bbox clipping) and prints the report.

## Running

Rendering needs matplotlib and numpy only — the `[plot]` extra of this project;
Figure S2 additionally uses pandas and numpy for the golden-file comparison.
The render path of S2 and S3 never imports `openreltime`, so the published
figures reproduce on a machine that does not have this package installed at
all.  Only `--regen-data` does (see below).

```bash
python figures/make_figures.py            --outdir /tmp/figs
python figures/附件图S1-软件架构图.py       --outdir /tmp/figs
python figures/附件图S2-性能对比图.py       --outdir /tmp/figs
python figures/附件图S3-置信区间图.py       --outdir /tmp/figs
```

**Always pass `--outdir`.**  The default output folder of every script is the
folder the script itself lives in, so a bare run writes its PDF/PNG (and
`.drawio` for Figure 1 and S1) into `figures/`.

Flags, per script:

| Script | `--outdir` | `--engine-root` | `--regen-data` | `--regen-out` | `--force` | other |
|---|---|---|---|---|---|---|
| `make_figures.py` | yes | — | — | — | — | — |
| `附件图S1-…​.py` | yes | — | — | — | — | `--figure-engine` |
| `附件图S2-…​.py` | yes | yes | yes | yes | yes | `--golden-dir` |
| `附件图S3-…​.py` | yes | yes | yes | yes | yes | `--golden-dir` (accepted, ignored) |

* **Figure S1 shares its drawing engine with Figure 1.**  It imports
  `make_figures.py`, looked up in `--figure-engine PATH`, in
  `$OPENRELTIME_FIGURE_ENGINE`, in its own folder, then in sibling folders and
  the working directory.  This folder holds both files, so no flag is needed.
* **`--engine-root`** names the OpenRelTime repository (the folder holding
  `openreltime/` and `data/`); it is also read from `$OPENRELTIME_ENGINE_ROOT`
  and probed by content, never by folder name.  When the scripts run from here
  the parent directory *is* the repository, so the flag is optional; it is
  needed only when the repository sits elsewhere.  `--golden-dir` similarly
  names `data/golden/r3f` (`$OPENRELTIME_GOLDEN_DIR`), which Figure S2 reads.

## Render environment of the published figures

The published PDF/PNG were rendered with **Python 3.11.15 and matplotlib
3.11.0** (numpy 2.4.6, pandas 2.3.3, Pillow 12.2.0) on macOS/arm64.  That is
verified rather than asserted: re-rendering all four scripts from this folder
reproduces the four published PNGs **byte for byte**, and the PDFs differ only
in their embedded creation timestamp.

Do not confuse this with the `environment` block inside
`附件图S2-性能对比图-data.json`, which records Python 3.14.6 / matplotlib 3.11.1
/ numpy 2.4.4 / pandas 3.0.3: those are the versions that *measured* the
wall-clock benchmark, not the versions that *drew* the figure.

## The checked-in data files are the published record

Figures S2 and S3 draw every number from the `-data.json` file checked in
beside them; the render path reads no other input.  `--regen-data` instead
**re-measures on the local machine** with the installed `openreltime`, and
because the checked-in file is the published record, a refresh writes a *copy*
next to it (`…-data.regen.json`) unless you name a file with `--regen-out PATH`
or consciously pass `--force`, which replaces it.

```bash
python figures/附件图S2-性能对比图.py --regen-data --regen-out /tmp/s2.json
python figures/附件图S3-置信区间图.py --regen-data --regen-out /tmp/s3.json
```

Read that guard as part of the data-availability statement.  `--regen-data`
is a legitimate way to re-check Figure S2 panel (a) — double-precision
agreement with R3F is machine-independent — and never a way to reproduce panels
(b)/(c), whose timings belong to the machine recorded in the file.

## Known limitation: `附件图S3-置信区间图.py --regen-data` aborts on 0.1.1

Re-measuring Figure S3 with the 0.1.1 engine **fails**.  Its panel (a) and
panel (c) inputs are `data/examples/example_calibrations.tsv`, whose three
placeholder rows are mutually inconsistent; 0.1.0 repaired that by rescaling
lineage rates and returned the result with warnings, whereas 0.1.1 checks the
converged ages against the tree and refuses them (`openreltime/calibrate.py`,
`_check_time_order`).  Measured on the current repository:

```
$ python 附件图S3-置信区间图.py --engine-root <repo> --regen-data --regen-out /tmp/x.json
...
ValueError: conflicting calibration constraints: the repair loop converged on
node ages that are not nested in the tree -- 1 node(s) older than the present
(negative age): [307]; 2 parent-child inversion(s): node 306 at 3225.03 sits
below its ancestor node 305 at 103.04; node 372 at 143.919 sits below its
ancestor node 371 at 6.23598.  This means the calibration set is unsatisfiable
by a rate rescaling of this tree; relax the bounds or move them toward the
root.
```

(The offending node ids and ages depend on the seed; the abort itself does not.)

`附件图S3-置信区间图-data.json` **predates that guard** — it records
`"openreltime_version": "0.1.0"` — and it *is* the published record: the
numbers in the accepted figure are exactly those in this file.  Nothing about
the figure needs to change, because panel (c) plots the very inconsistency that
the check now rejects; only the refresh path is affected.  `data/examples/README.md`
explains why the file is a format illustration that is expected to fail, and
`data/examples/example_calibrations_feasible.tsv` is the runnable counterpart.

## Distribution

`MANIFEST.in` ships `figures/*.py`, `figures/*.json` and this README with the
**source distribution**, on the same pattern it already uses for `data/`,
`docs/` and `reproduce/`.  Those directories are not Python packages, so
setuptools cannot and does not place them inside the wheel: the installed
layout is unchanged (`openreltime` plus the `openreltime`/`ort` console
scripts), and the figure sources reach a reader through the sdist or a
checkout.  Rendered artwork is excluded, so a local `--outdir figures` run
cannot leak into a distribution built from the working tree.

## Style

The four scripts and two data files are **byte-identical copies** of the
sources that produced the published render, so they are not reformatted to this
repository's ruff profile; `[tool.ruff].extend-exclude` lists `figures` and the
lint gate therefore stays at its pre-existing count.  They are also outside the
package import graph and outside the test suite.
