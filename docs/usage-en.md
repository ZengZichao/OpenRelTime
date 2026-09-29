# OpenRelTime user guide

**Language: English** · [中文](usage-zh.md)

Version 0.1.0 · Requires Python ≥ 3.10 · Licence GPL-3.0-or-later

This guide covers installation, input preparation, every subcommand, the
calibration file format, the output inventory, a pointer to the separately
installed desktop GUI, FAQs and troubleshooting.  For the Chinese version see
[usage-zh.md](usage-zh.md); for the semantics of
every option see [parameters.md](parameters.md); for algorithms and
limitations see [methods.md](methods.md).

---

## 1. Installation

### 1.1 From a source checkout

> **Distribution.** This release is distributed as source. Install with
> `pip install --editable ".[dev,plot]"` from a checkout; an index release will
> be announced here.

```bash
git clone https://github.com/ZengZichao/OpenRelTime
cd OpenRelTime
python -m pip install --editable ".[dev,plot]"   # core + plotting + dev tools
```

`plot` (matplotlib visualisation) and `dev` (pytest, pytest-cov, ruff, mypy for
contributors) are the only extras the package declares
(`[project.optional-dependencies]` in `pyproject.toml`); the
analytical confidence intervals are pure SciPy/NumPy and need no extra.  Combine
extras in one call when you need several, e.g.
`pip install --editable ".[plot,dev]"`.  The
desktop GUI is **not** an extra of this package — it is its own distribution,
`OpenRelTime-Studio`, which depends on this engine (§7).
Verify with `openreltime --version` (short alias `ort`).  In a
conda/micromamba environment, activate the target environment first so the
package lands in it:

```bash
micromamba activate <your-env>      # or: conda activate <your-env>
pip install --editable ".[dev,plot]"
```

There are no
compiled extensions; installation is fast on Linux, macOS and Windows.

### 1.2 Working against a checkout

An **editable install is the supported way** to run analysis, benchmark or
figure scripts that live outside the repository against a source checkout: the
checkout (and therefore the `openreltime` package and the console scripts
`openreltime`/`ort`) is then importable from any working
directory, so a script simply writes

```python
import openreltime as ort
```

Do **not** prepend the checkout to `sys.path` by hand (e.g.
`sys.path.insert(0, "/path/to/OpenRelTime-source")`) — that hack bypasses the
declared dependency set, breaks the console scripts, and silently shadows the
installed package.  Because the install is editable, edits to the checkout take
effect immediately without reinstalling.  Confirm what a script actually
resolves to with:

```bash
python -c "import openreltime; print(openreltime.__version__, openreltime.__file__)"
```

## 2. Input preparation

### 2.1 Phylogeny (required)

A **rooted binary** Newick (or NEXUS) tree whose branch lengths are expected
substitutions per site — the `.treefile` of IQ-TREE/RAxML/FastTree works
directly.  Requirements:

- at least three ingroup taxa;
- finite, non-negative branch lengths (missing lengths become 0 with a
  warning);
- polytomies are refused by default; `--resolve random` binarises them with
  zero-length branches (logged).

### 2.2 Outgroup (recommended)

The tree is rooted on the outgroup MRCA and the outgroup tips are removed;
all outputs contain ingroup only.  Two forms are supported:

- **single-tip outgroup**: one terminal name;
- **multi-tip outgroup**: several names (comma-separated or one per line in
  a text file).

```bash
--outgroup Ornithorhynchus_anatinus
--outgroup "Ornithorhynchus_anatinus,Zaglossus_bruijni,Tachyglossus_aculeatus"
--outgroup outgroup.txt      # one name per line
```

**Omitting `--outgroup`.**  The flag is optional, and then no rooting and no
pruning happen: the tree is used exactly as read.  Every RelTime quantity is
anchored on the root, so the root position in the input file silently determines
all node ages and OpenRelTime does **not** warn about it — confirm the input is
rooted on the clade you intend (R3F requires an outgroup; see
`docs/validation.md` §5).

**Complete-clade check.**  A multi-tip outgroup must form a complete clade
(a monophyletic group): after tracing the specified branches back to their
most recent common ancestor (MRCA), *every* tip descending from that ancestor
must belong to the outgroup.  If the MRCA also descends to unspecified tips,
the run fails with an error listing them; `--outgroup-check warn` downgrades
this to a warning (rooting still uses the MRCA).  A single-tip outgroup is
always complete and is not checked.

**Outgroup by group name.**  Each `--outgroup` token is first matched
against tip names; otherwise it is resolved as a group name (§2.4), expanded
to its member tips and subjected to the same complete-clade check:

```bash
--outgroup Monotremata --taxon-table taxa.tsv
```

### 2.3 Calibrations (for `calibrate`)

A tab-separated file (§4) or a list of `openreltime.Calibration` objects.
`taxon_set` entries may be tip labels or a group name (resolved via
`--taxon-table` or auto-detected label prefixes, §2.4); a group judged
non-monophyletic raises a warning and its MRCA is used as the calibration
node.

### 2.4 Tip classification (monophyly assessment)

The `monophyly` subcommand, group-name outgroups and group-name calibrations
need a tip → group mapping, taken from either source:

1. **user-supplied table** (`--taxon-table`): TSV/CSV with two columns; the
   headers may be any of `tip`/`label`/`name` and `group`/`clade`/`taxa`/
   `classification` (example: `data/examples/example_taxa.tsv`);
2. **automatic recognition from tip labels** (default): following the
   ``Genus_species`` convention, the text before the first separator
   (`_`, space, `|`, `@`) is taken as the group (genus), e.g.
   `Mus_musculus` → group `Mus`.

## 3. Command-line usage

### 3.1 Relative rates and times

```bash
openreltime rates       -i tree.nwk --outgroup og.txt -o mam
openreltime times       -i tree.nwk --outgroup og.txt -o mam --r3f-compat
openreltime rates-times -i tree.nwk --outgroup og.txt -o mam --plot mam_tree.png
```

- `--normalize` scales the maximum node time to 1.0 (MEGA/PNAS convention);
- `--mean geometric|arithmetic` selects the averaging convention; geometric is
  R3F-compatible and the default everywhere;
- `--no-guard` disables the extreme-rate guard (default threshold 20, as in
  R3F), `--rate-ratio-threshold FLOAT` sets it.  The guard belongs to the
  **time** pass only: `rates` accepts `--rate-ratio-threshold` for signature
  symmetry, ignores it, and records a warning saying so;
- `--r3f-compat` writes R3F file names and column names
  (`_RRF_times.csv`, …) so existing pipelines can swap tools.

`--fmt`, `--resolve`, `--taxon-table` and `--outgroup-check` (§2) are available
on the tree-based subcommands.  `docs/parameters.md` opens with a matrix of
which flag each subcommand takes.

### 3.2 Calibration to absolute times

```bash
openreltime calibrate -i tree.nwk -c calibrations.tsv --outgroup og.txt \
    --method effective --n-effective 10000 --seed 42 -o mam_cal
```

- `--method bounds` (default): hard min/max constraints.  The global factor
  `f` is the midpoint of the feasible intersection; residual violations are
  fixed by proportional rate scaling along the offending lineage
  (propagated to descendants, ≤ 100 iterations), otherwise the conflict is
  reported;
- `--method effective`: msz236 density resampling; `--n-effective` sets the
  replicate count (default 10,000).  Outputs include the effective bounds
  table and empirical node-age quantiles.  A density-only row needs this
  method — under `bounds` it constrains nothing and is rejected at parse time;
- `--mean`, `--rate-ratio-threshold` and `--no-guard` choose the RRF estimate
  being calibrated; whatever you pass is recorded in the report and replayed by
  `openreltime ci`, so the intervals belong to the same estimator.

### 3.3 Analytical confidence intervals

```bash
openreltime ci -c mam_cal --n-sites 1000 -o mam_ci
```

- `-c` is the **output prefix** of the `calibrate` run, not a tree: the tree,
  outgroup and input handling are re-read from `<prefix>_report.json`;
- `--branch-var var.tsv` — user-provided sampling variance per edge
  (exactly the columns `node_id`, `var`);
- `--n-sites L` — Poisson approximation `vS(b) = b/L`;
- with neither, the interval captures rate-heterogeneity variance only (the
  msz236 simulation protocol);
- `--level` sets the confidence level; intervals are truncated at hard
  bounds and clipped at 0;
- `--mean`, `--rate-ratio-threshold` and `--no-guard` override the settings
  `calibrate` recorded.  Left unset, the intervals are differentiated through
  the same estimator that produced the point estimate — including
  `--method effective` results, which are differentiated along their effective
  bounds;
- `<prefix>_ci.csv` lists `node_id, label, time, se, lower, upper, width,
  se_reliable, notes`.  A node whose age is pinned by its own calibration has a
  zero delta-method derivative, so it is reported with the uncertainty that
  calibration asserts and flagged `se_reliable=False`; `notes` says why, and
  `n_nodes_se_not_reliable` counts them in the report JSON.

### 3.4 CorrTest

```bash
openreltime corrtest -i tree.nwk --outgroup og.txt -o mam
openreltime corrtest -i tree.nwk --outgroup og.txt --sister-resample 100 --seed 1 -o mam
```

`<prefix>_corrtest.txt` reports CorrScore ∈ [0, 1] and the P band
(≥ 0.92 → P < 0.001; [0.83, 0.92) → P < 0.01; [0.5, 0.83) → P < 0.05;
< 0.5 → P > 0.05).  Scores ≥ 0.5 indicate significant rate autocorrelation.

### 3.5 ddBD tree prior

```bash
openreltime ddbd -i tree.nwk --outgroup og.txt --anchor-time 1.85 -o mam_bd
openreltime ddbd -i tree.nwk --sampling-frac 0.5 -o mam_bd
```

- `--anchor-node N --anchor-time T` calibrates one internal node age (keep
  the maximum age ≤ 10 time units);
- `--measure SSE|KL` selects the initial-value grid score;
- `--sampling-frac` sets the sampling fraction to **report**: as in R3F the
  likelihood fit always leaves rho free, so the value replaces the reported
  fraction without constraining the optimisation (the free-rho optimum stays
  available as `sampling_frac_fitted` in the report);
- `--plot curve.png` saves the node-time histogram with the fitted
  birth-death density.

### 3.6 tree2table

```bash
openreltime tree2table -i tree.nwk -o mam_tbl          # branch lengths
openreltime tree2table -i tree.nwk --time -o mam_tbl   # node ages
```

### 3.7 Monophyly assessment

Test whether taxon groups are monophyletic (complete clades).  Groups come
from `--taxon-table` or are auto-detected from tip-label prefixes; `--groups`
selects which groups to test (default: all) and `--tips` tests an arbitrary
tip set:

```bash
# auto-detect genus groups from labels, test all
openreltime monophyly -i tree.nwk -o chk

# test named groups from a classification table
openreltime monophyly -i tree.nwk --taxon-table taxa.tsv --groups Monotremata,Primates -o chk

# test an arbitrary tip set
openreltime monophyly -i tree.nwk --tips A,B,C -o chk
```

Output `<prefix>_monophyly.csv`: one row per group with tip counts, the MRCA
tip count, extra tips under the MRCA and the `is_monophyletic` verdict.
Non-monophyletic groups additionally raise a warning on stderr.  The same
verdict backs group-name outgroups (§2.2) and group-name calibrations (§4).
Without `--taxon-table`, groups are auto-detected from tip-label prefixes;
`--group-sep REGEX` changes the split pattern (default `[_|@\s]`, i.e. the
first of `_`, `|`, `@` or whitespace).

### 3.8 BLB pipeline (needs IQ-TREE)

```bash
openreltime blb -a alignment.fasta --iqtree /path/to/iqtree \
    --outgroup og.txt --gamma 0.7 --n-subsample 10 --n-replicate 10 \
    --iqtree-arg=-m --iqtree-arg=GTR+G --seed 42 -o blb_out
```

Two-level resampling: (i) `ceil(L**gamma)` sites without replacement;
(ii) with-replacement resampling back to length L; IQ-TREE infers each
replicate tree and OpenRelTime times it.  `blb_out/blb_summary.csv` reports
per-clade medians, 2.5/97.5 percentiles and the `topology_ok` fraction, and
`blb_replicates.csv` keeps every replicate table.  `--workdir` retains the
intermediate replicate alignments and trees (by default they go to a temporary
directory that is removed); `--seed` also reaches each IQ-TREE run as `-seed`.
The
`--gamma 0.7` / `10 x 10` defaults sit inside the admissible BLB range but do
**not** reproduce the published RelTime-JA run (`g = 0.78`, `20 x 20`); see
`docs/parameters.md` and `docs/validation.md`.

### 3.9 MEGA-CC bridge (optional)

```bash
openreltime megacc -i tree.nwk -c calibrations.tsv --megacc /path/to/megacc \
    --outgroup og.txt --workdir megacc_run -o mam_mega
```

Generates `reltimeFromBranchLengths.mao`, runs megacc and parses its output
back into `<prefix>_megacc_times.csv`, an R3F-compatible node-time table in
which a `source` column records where every row came from.  `-c/--calibrations`
is optional (branch-lengths-only RelTime without it) and `--mao custom.mao`
replaces the generated control file.  If megacc produces nothing parseable the
run raises — it never substitutes OpenRelTime's own estimates and labels them
as MEGA's.
Two limits apply: the bridge needs a real `megacc` executable, which CI neither
installs nor runs, so MEGA-CC agreement is established against the published
R3F/MEGA results and the shipped golden files rather than by a head-to-head run
here; and
OpenRelTime's analytical RRF layer lacks MEGA's likelihood-ratio refinement, so
node-by-node differences from MEGA are expected and this is a comparison tool,
not an equivalence claim (see `docs/validation.md`).

## 4. Calibration file format

Tab-separated (example: `data/examples/example_calibrations.tsv` — a **format
illustration only**; the numbers are arbitrary illustrative values, see
`data/examples/README.md`).  Example and
golden data ship with the source distribution and an editable checkout
(`data/`, `MANIFEST.in`), not with a pure wheel: run the examples from the
checkout.  Header row required; every row carries all six tab-separated
fields (`.` marks an empty cell).

| Column | Required | Meaning |
|---|---|---|
| `node_id` | one of the two | internal node id (as in `_times.csv`) |
| `taxon_set` | one of the two | `\|`-separated tip labels, or a single group name (`--taxon-table` / auto-detection, §2.4); target = their MRCA; a non-monophyletic group raises a warning |
| `min_bound` | ≥ 1 bound or a density | minimum age; `.` = missing |
| `max_bound` | | maximum age |
| `density` | optional | `uniform` / `exponential` / `normal` / `lognormal` |
| `density_params` | with density | `key=value;…`, e.g. `offset=60;mean=20` |

Parameter keys — uniform `min,max`; exponential `offset,mean`; normal
`mean,sd`; lognormal `offset,meanlog,sdlog`.  Values are parsed strictly;
no expression evaluation.

## 5. Output inventory

| Suffix | Written by | Content |
|---|---|---|
| `_rates.csv` | `rates` | NodeLabel, NodeId, Des1, Des2, Rate |
| `_rates.nwk` | `rates` | same tree with rates as branch annotations |
| `_times.csv` | `times` | …, Time |
| `_timetree.nwk` | `times` | relative-time tree |
| `_rates_times.csv` | `rates-times` | Rate + Time |
| `_timetree.nexus` | `rates-times` | with `[&rate=...]`, opens in FigTree |
| `_RRF_*` variants | `--r3f-compat` (`times`, `rates-times`) or the API `write(r3f_compat=True)` | `_RRF_rates.csv`, `_RRF_times.csv`, `_RRF_timetree.nwk`, `_RRF_table.csv`, `_RRF_timetree.nexus` |
| `<prefix>.csv` | `tree2table` | node table (branch lengths, or ages with `--time`) |
| `_monophyly.csv` | `monophyly` | group, n_tips, is_monophyletic, n_mrca_tips, extra_tips, missing_names, tips |
| `_calibrated.csv/.nwk/.nexus` | `calibrate` | absolute-time table and trees |
| `_effective_bounds.csv` | `calibrate --method effective` | effective bounds (msz236) + empirical quantiles |
| `_ci.csv` | `ci` | node_id, label, time, se, lower, upper, width, se_reliable, notes |
| `_corrtest.txt` / `_ddbd.txt` | `corrtest` / `ddbd` | CorrScore / BD parameters |
| `blb_summary.csv`, `blb_replicates.csv` | `blb` | per-clade summary; every replicate table |
| `_megacc_times.csv` | `megacc` | MEGA-CC node times with a `source` column |
| `_report.json` | see below | all parameters, warnings and the random seed |

`_report.json` is written by the analysis stages whose parameters can change
the numbers: `rates`, `times`, `rates-times`, `calibrate`, `ci`, `corrtest`
and `ddbd`.  Runs that share a prefix **merge** into one report file, each
contributing its own top-level key, so a whole pipeline stays inspectable in a
single JSON.  `tree2table`, `monophyly`, `blb` and `megacc` write only their
table.

## 6. Python API at a glance

```python
import openreltime as ort

tree  = ort.read_tree(path, fmt="newick", outgroup=[...])
rates = ort.rrf_rates(tree, mean="geometric", rate_ratio_threshold=None)
times = ort.rrf_times(tree, rate_ratio_threshold=20.0, normalize=False)
rt    = ort.rrf_rates_times(tree)
cal   = ort.calibrate(times, cals, method="effective", n_effective=10_000, seed=42)
ci    = ort.confidence_interval(cal, seq_length=1000, level=0.95)
corr  = ort.corrtest(tree, sister_resample=100, seed=1)
bd    = ort.ddbd(tree, anchor_time=1.85, measure="SSE")
table = ort.tree2table(tree, time=False)

# monophyly assessment (classification table or auto-detected label prefixes)
table_map = ort.parse_taxon_table("taxa.tsv")                # {group: [tips]}
verdict = ort.check_monophyly(tree, ["Mus_musculus", "Mus_spretus"], name="Mus")
verdict.is_monophyletic
groups  = ort.detect_groups_from_labels(tree.tip_labels())   # {genus: [tips]}
# group-name outgroup / calibration:
tree2 = ort.read_tree(path, outgroup=["Monotremata"], taxon_table=ort.parse_taxon_table("taxa.tsv"))
cal2  = ort.calibrate(times, cals, taxon_table=table_map)
```

All result objects expose `to_pandas()`, `to_json()` and `write(prefix)`;
API and CLI outputs are produced by the same code path.

## 7. Desktop GUI — OpenRelTime Studio (separate product)

The point-and-click front end for this engine, **OpenRelTime Studio**, is not
part of this repository and is not installed from it.  It ships as its own
distribution and depends on this engine's public API, exactly as any third-party
script does — so every number it produces is the number you would get from the
CLI or the API with the same inputs (the analysis code itself lives here, and
nothing about it changes when you add or remove the GUI).

```bash
pip install OpenRelTime-Studio     # pulls in this engine as a dependency
openreltime-studio
```

* Home page: <https://github.com/ZengZichao/OpenRelTime-Studio>
* Interface capabilities: a fully bilingual UI (English / 简体中文, switched live from
  *View ▸ Language*), light and dark themes (*View ▸ Theme*), four root orientations
  (*View ▸ Root position*), a hand-drawn vector SVG icon, and a built-in example tree
  with demo calibrations (*File ▸ Open Bundled Example…*).
* What it does: load a Newick/NEXUS tree, run RRF, add calibrations by clicking
  internal nodes on the canvas, then run calibration, analytical CIs, CorrTest
  and ddBD; export CSV/NEXUS/JSON/PNG plus the equivalent command line.
* Window layout, the run-button prerequisite chain, what export writes and the
  rest of the GUI manual: **see that project's own documentation**.  This guide
  documents the engine only.
* The calibration file it loads is the same format as §4, and the files it
  exports are the same as §5.

## 8. FAQ

**Q1. Are results identical to MEGA's RelTime?**
No.  OpenRelTime implements the analytical RRF layer (the layer R3F also
implements); MEGA adds an LRT-based rate-constraint refinement.  Expect
correlations of ~0.98–0.99, not bitwise equality — see `docs/methods.md`.

**Q2. How close is it to R3F?**
The G2a/G2b gate asserts a regression slope >= 0.999 and a median relative
deviation <= 1e-6 on the 274-tip mammalian benchmark; the observed medians are
~6.8e-16 (rates) and ~8.7e-16 (times), with maximum deviations ~4.8e-15 and
3.6e-14 that are measured but not gated (see `docs/validation.md`).

**Q3. The tree has many zero-length branches.**
OpenRelTime warns when ≥ 10% of branches are zero (the R3F threshold).
Affected estimates are unreliable; consider re-inferring the tree.  The
epsilon handling matches R3F exactly, so sensitivity comparisons are
possible.

**Q4. Polytomies raise an error.**
By design — the analytical solution assumes binary trees.  Use
`--resolve random` after checking the topology; affected nodes are logged.

**Q5. BLB is slow.**
BLB runs `n_subsample × n_replicate` IQ-TREE jobs; the OpenRelTime timing
step itself is sub-second per tree (1000 tips < 1 s).  Budget accordingly.

**Q6. How to cite?**
See `CITATION.cff`; please also cite msy044 (RRF), msz236 (CIs/effective
bounds), msz014 (CorrTest) and btab307 (ddBD).

## 9. Troubleshooting

| Message / symptom | Cause and fix |
|---|---|
| `tree contains polytomies at ...` | binarise first or `--resolve random` |
| `negative branch length at node ...` | fix the tree file; negative lengths are invalid |
| `duplicate tip names` | rename duplicated taxa |
| `ingroup has N tips after outgroup removal` | fewer than 3 ingroup taxa remain; adjust the outgroup |
| `outgroup [...] is not a complete clade` | the multi-tip outgroup MRCA also descends to unspecified tips; complete the outgroup or use `--outgroup-check warn` |
| `outgroup taxon [...] not found in tree (as tip or group)` | the name is neither a tip nor a group; check the spelling or pass `--taxon-table` |
| `conflicting calibration constraints` | calibration nodes contradict each other (e.g. child bounds older than parent); relax bounds |
| `could not be satisfied within 100 iterations` | same; see the listed nodes |
| `a calibration needs at least one finite bound … but neither was given` | the row has `.` for `min_bound` and `max_bound`; add hard bounds, or run `--method effective` for a density-only row |
| `Under method='bounds' a density does not count as a bound` | density-only calibration with `--method bounds`; switch to `--method effective` |
| `no 'calibrate' block in <prefix>_report.json` | `ci -c` was given a prefix that never ran `calibrate`; re-run `calibrate` with that `-o` |
| `cannot locate the original tree` | the tree file moved since `calibrate` ran; re-run `calibrate` with the same `-i`/`-o` so `tree_file` is recorded |
| `group '…' not found in …; available: …` | the name is neither in `--taxon-table` nor an auto-detected prefix; the message lists what is available |
| `sampling_frac must lie in (0, 1]` | `0` is R3F's "estimate it" sentinel — omit `--sampling-frac` instead |
| `anchor node id must be a positive internal-node id` | node ids are 1-based here; `0` is not the "no anchor" sentinel |
| `executable not found: iqtree/megacc` | pass the full executable path |
| `argument ... contains characters outside the allowed set` | external-tool arguments must not contain shell metacharacters or spaces |
| CorrScore = 0 or NaN | tree too small or extreme rates; try `--sister-resample`, then inspect branch lengths |

## 10. Reproducibility

- Golden standards and the generating script: `data/golden/r3f/`,
  `reproduce/golden_r3f.R` (R 4.5.3, ape 5.8.1, phangorn 2.11.1, R.utils
  2.13.0, FNN 1.1.4.1, RColorBrewer 1.1-3; see `data/PROVENANCE.md`);
- tests: `python -m pytest` from the repository root runs the engine suite
  (`tests/`).  The desktop front end carries its own tests in the separate
  `OpenRelTime-Studio` project (§7);
- continuous integration: `.github/workflows/ci.yml` runs that suite on
  Python 3.10, 3.11, 3.12 and 3.13 on Ubuntu and macOS.  Windows is a declared
  target but is **not** covered by CI, and external-tool paths (IQ-TREE for
  `blb`, MEGA-CC for `megacc`) are not installed in the workflow — see
  `docs/validation.md` for the tested-vs-claimed breakdown;
- all randomness (effective bounds, BLB, sister resampling) flows through
  `--seed` and is recorded in `_report.json`.
