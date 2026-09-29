# Parameter handbook

**Language: English** · [中文](parameters-zh.md)

Every optional parameter of OpenRelTime: semantics, default, valid range and
recommendations.

Legend — **API**: keyword argument; **CLI**: flag; **where**: module.

---

## Which flags each subcommand accepts

Read this table from `openreltime <subcommand> --help`; it is the quick index
into the entries below.  `-i/--input`, `-o/--output` and `-h/--help` exist on
every subcommand except `blb` (which takes `-o/--output-dir`) and are not
repeated.

| Flag | rates | times | rates-times | calibrate | ci | corrtest | ddbd | tree2table | monophyly | blb | megacc |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `--outgroup` | ✓ | ✓ | ✓ | ✓ | –¹ | ✓ | ✓ | – | – | ✓ | ✓ |
| `--fmt` | ✓ | ✓ | ✓ | ✓ | –¹ | ✓ | ✓ | ✓ | ✓ | – | – |
| `--resolve` | ✓ | ✓ | ✓ | ✓ | –¹ | ✓ | ✓ | – | ✓ | – | – |
| `--taxon-table` | ✓ | ✓ | ✓ | ✓ | –¹ | ✓ | ✓ | – | ✓ | – | – |
| `--outgroup-check` | ✓ | ✓ | ✓ | ✓ | –¹ | ✓ | ✓ | – | – | – | – |
| `--mean` | ✓ | ✓ | ✓ | ✓ | ✓² | – | – | – | – | – | – |
| `--rate-ratio-threshold` | ✓³ | ✓ | ✓ | ✓ | ✓² | – | – | – | – | – | – |
| `--no-guard` | – | ✓ | – | ✓ | ✓² | – | – | – | – | – | – |
| `--normalize` | – | ✓ | ✓ | – | – | – | – | – | – | – | – |
| `--r3f-compat` | – | ✓ | ✓ | – | – | – | – | – | – | – | – |
| `--plot` | – | – | ✓ | – | – | – | ✓ | – | – | – | – |
| `--calibrations` | – | – | – | ✓ | – | – | – | – | – | – | ✓⁴ |
| `--method` | – | – | – | ✓ | – | – | – | – | – | – | – |
| `--n-effective` | – | – | – | ✓ | – | – | – | – | – | – | – |
| `--seed` | – | – | – | ✓ | – | ✓ | – | – | – | ✓ | – |
| `--branch-var` | – | – | – | – | ✓ | – | – | – | – | – | – |
| `--n-sites` | – | – | – | – | ✓ | – | – | – | – | – | – |
| `--level` | – | – | – | – | ✓ | – | – | – | – | – | – |
| `--sister-resample` | – | – | – | – | – | ✓ | – | – | – | – | – |
| `--anchor-node` / `--anchor-time` | – | – | – | – | – | ✓ | ✓ | – | – | – | – |
| `--sampling-frac` | – | – | – | – | – | – | ✓ | – | – | – | – |
| `--measure` | – | – | – | – | – | – | ✓ | – | – | – | – |
| `--time` | – | – | – | – | – | – | – | ✓ | – | – | – |
| `--groups` / `--tips` / `--group-sep` | – | – | – | – | – | – | – | – | ✓ | – | – |
| `--alignment` / `--iqtree` / `--gamma` / `--n-subsample` / `--n-replicate` / `--iqtree-arg` | – | – | – | – | – | – | – | – | – | ✓ | – |
| `--megacc` / `--mao` | – | – | – | – | – | – | – | – | – | – | ✓ |
| `--workdir` | – | – | – | – | – | – | – | – | – | ✓ | ✓ |

¹ `ci` does not take tree inputs on the command line: it re-reads the tree,
outgroup, format, polytomy handling and outgroup check recorded by the
`calibrate` run in `<prefix>_report.json` (see `calibrated`).
² overrides the value recorded by `calibrate` (see `mean` /
`rate_ratio_threshold`); unset means "use what was recorded".
³ accepted but **ignored** by the rate pass — see `rate_ratio_threshold`.
⁴ optional for `megacc`: without it the bridge runs RelTime on branch lengths
alone.

---

## Input handling

### `fmt` — input format
* Values: `"newick"` | `"nexus"`.
* Default: `"newick"` (`openreltime.read_tree(path, fmt="newick")`; the CLI
  flag `--fmt` also defaults to `newick` and accepts only these two values).
* NEXUS files may carry `[&rate=...]` annotations (R3F/FigTree dialect);
  they survive reading instead of being stripped, and are regenerated on
  write (gate G2c asserts the round trip on the R3F golden NEXUS).

### `outgroup` — rooting outgroup
* Type: list of tip names (API), comma-separated string or one-per-line
  text file (CLI: `--outgroup`).
* Default: `None` — no rooting and no pruning; the tree is used exactly as
  read and must already be rooted and binary.  The CLI flag is likewise
  optional and unset means `None`.
* The tree is rooted on the outgroup MRCA and the outgroup tips are removed;
  output contains ingroup tips only.  Behaviour matches
  `ape::root(resolve.root=TRUE)` + `ape::drop.tip` (ADR-004).
* Recommendation: use 1-5 slowly evolving taxa that form a clade.
* Two input forms: a **single tip**, or **multiple tips** that must form a
  complete clade (see `outgroup_check`).  A token that is not a tip name is
  resolved as a **group name** via `taxon_table` or auto-detected tip-label
  prefixes and expanded to its members.

### `taxon_table` — tip classification
* Type: `None` | `{group: [tip, ...]}` mapping (API:
  `openreltime.parse_taxon_table(path)`; CLI: `--taxon-table PATH`).
* TSV/CSV with a tip column (`tip`/`label`/`name`/…) and a group column
  (`group`/`clade`/`taxa`/`classification`/…); with two unrecognised
  headers, first column = tip, second = group.  Without a table, groups are
  auto-detected from tip-label prefixes (`Genus_species` → `Genus`).
* Used by group-name outgroups, group-name calibration `taxon_set` values
  and the `monophyly` subcommand.

### `outgroup_check` — multi-tip outgroup completeness
* Values: `"error"` (default) | `"warn"`.  CLI: `--outgroup-check`.
* A multi-tip outgroup must be a complete clade: all tips under the outgroup
  MRCA must belong to the outgroup.  `"error"` raises `ValueError` listing
  the extra tips; `"warn"` logs a warning and roots on the MRCA anyway.
  Single-tip outgroups are always complete and are never checked.

### `monophyly` subcommand
* `openreltime monophyly -i tree.nwk [--taxon-table taxa.tsv] [--groups A,B]
  [--tips A,B,C] [--group-sep REGEX] [--fmt ...] [--resolve ...] -o PREFIX`.
* Tests whether each group is monophyletic — i.e. the MRCA of its tips
  descends to no other tip (the same criterion as `outgroup_check`).
* `--groups` selects which groups to test (default: every group in the
  classification); `--tips` tests one ad-hoc tip set instead; the two are
  independent ways of narrowing the run.
* `--group-sep REGEX` — separator pattern for the automatic label-prefix
  recognition; default splits on the first of `_`, `|`, `@` or whitespace.
  Only relevant without `--taxon-table`.
* Writes `<prefix>_monophyly.csv` (columns `group`, `n_tips`,
  `is_monophyletic`, `n_mrca_tips`, `extra_tips`, `missing_names`, `tips`);
  non-monophyletic groups log a warning on stderr.  `missing_names` lists
  requested names absent from the tree.

### `resolve_polytomy` — non-binary trees
* Values: `"error"` (default) | `"random"`.  CLI: `--resolve error|random`.
* `"random"` resolves polytomies with zero-length branches (a warning is
  logged).  Zero branches make RRF rates/times on affected nodes unreliable
  (R3F warns above 10% zero-length branches; OpenRelTime logs the same
  threshold).

## RRF estimation

### `mean` — averaging convention
* Values: `"geometric"` (default) | `"arithmetic"`.
* Geometric = msy044 eqs 28-42, R3F-compatible.  Arithmetic = PNAS 2012
  semantics (eqs 1-27), provided for cross-validation; outputs will differ.

### `rate_ratio_threshold` — extreme-rate guard
* Type: float | `None`.  Default: `None` in `rrf_rates`, `20.0` in
  `rrf_times`/`rrf_rates_times`.  CLI: `--rate-ratio-threshold` /
  `--no-guard` (see the note on `rates` below).
* After the preorder adjustment, a node whose adjusted rate exceeds the
  threshold (or falls below its reciprocal) has its node time replaced by
  the closest non-exceeding ancestor's time (R3F semantics).  `None`
  disables the guard.
* **The guard belongs to the time pass only.**  `rrf_rates` takes the
  keyword for signature symmetry with `rrf_times` but never applies it: the
  R3F reference guards only where it rewrites node times (`rrf_times.R`
  lines 353-373), and `rrf_rates.R` has no times to rewrite.  Passing a
  number to `rrf_rates` therefore logs a warning, which is also recorded in
  `RateResult.warnings` and in `<prefix>_report.json`.  `openreltime rates
  --rate-ratio-threshold` reaches that code path unchanged.
* Sensitivity: repeating an analysis with `None`/20/10 quantifies how many
  nodes are affected (`n_rate_guarded_nodes` in the report JSON).
* `openreltime ci` re-runs the estimator behind the intervals, so it exposes
  `--mean`, `--rate-ratio-threshold` and `--no-guard` to override what the
  `calibrate` run recorded; leaving them unset reproduces the recorded
  settings exactly.  `--no-guard` and `--rate-ratio-threshold` are mutually
  exclusive.

### `normalize` — rescale node times
* Default: `False`.  CLI: `--normalize`.
* Divides node times by their maximum (root time = 1.0), the MEGA/PNAS 2012
  reporting convention.  R3F reports raw adjusted times (default here too).

## Calibration (`calibrate`)

### `method`
* `"bounds"` (default) — use min/max constraints directly.  The global time
  factor `f` is the midpoint of the intersection of the feasible intervals
  `[min/t, max/t]`; residual violations trigger proportional rate scaling of
  the offending lineage, propagated to descendants, up to 100 iterations
  (`openreltime._constants.CALIBRATION_MAX_ITER`).  Infeasible constraint
  sets raise `ValueError` with per-node diagnostics, and a feasible factor
  that comes out non-finite raises rather than propagating `f = inf` into
  every absolute age.  A density-only row constrains nothing under
  `method="bounds"`, so it is rejected when the calibration table is parsed
  (`--method effective`, or add `min_bound`/`max_bound`).
* `"effective"` — msz236: for each density calibration two dates are sampled
  per replicate as (min, max); the analysis is repeated `n_effective` times;
  the 2.5/97.5 percentiles per calibrated node become *effective bounds*
  used for the final estimate.  Empirical 2.5/97.5 node-age quantiles are
  reported alongside.

### `n_effective`
* Default `10_000` (msz236 protocol); lower it (>= 200) for exploratory runs
  on large trees.  Replicates reuse the cached composite lengths, so each
  replicate costs O(#calibrations + n).

### `seed`
* Randomness (density sampling, BLB, sister resampling) is always routed
  through a `numpy.random.Generator`; record the seed (stored in
  `<prefix>_report.json`) for reproducibility.

### Calibration file format (`calibrations.tsv`, TSV)
Columns: `node_id` (optional), `taxon_set` (`|`-separated, optional;
exactly one of the two), `min_bound`, `max_bound` (`.` = missing),
`density` (optional: `uniform|exponential|normal|lognormal`),
`density_params` (`key=value` pairs separated by `;`):
* uniform: `min,max`
* exponential: `offset,mean`
* normal: `mean,sd`
* lognormal: `offset,meanlog,sdlog`

`taxon_set` entries may be tip labels or a single **group name** (resolved
against `taxon_table` or auto-detected tip-label prefixes).  Whenever a
`taxon_set` turns out not to be monophyletic on the ingroup tree, a warning
is logged and the MRCA of its tips is used as the calibration node.

Parameters are parsed strictly as key/value numbers; expression evaluation
is never performed.

**Every row must constrain something.**  A row carrying neither
`min_bound` nor `max_bound` (and no density) is rejected at parse time with
the offending line number, because it would otherwise drive the global factor
`f` to `inf` and turn every absolute age into `inf` silently.  A density-only
row is likewise rejected under `method="bounds"`; use `method="effective"`
for it, or supply hard bounds.  Missing columns, unparseable numbers and
duplicate node ids raise too.

## Confidence intervals (`confidence_interval`)

### `calibrated` / CLI `-c`, `--calibrated`
* The input is a `CalibratedResult` (API) or the **output prefix of a
  previous `openreltime calibrate` run** (CLI).  The CLI re-reads the tree,
  outgroup, input format, polytomy handling and outgroup check recorded in
  `<prefix>_report.json`, so the intervals belong to the estimate being
  qualified; it fails with a readable message when the report or the original
  tree file cannot be found.
* One estimator throughout: the differentiation pipeline reproduces the
  averaging convention (`mean`), the guard threshold
  (`rate_ratio_threshold`) and the calibration bounds **actually used** for
  the point estimate.  A `method="effective"` result is differentiated along
  its effective bounds (msz236 treats them as the constraints), a
  `mean="arithmetic"` result along arithmetic derivatives.  A `method` other
  than `bounds`/`effective` raises, and an effective result that carries no
  effective bounds (e.g. `n_effective < 2`) raises too.
* Output `<prefix>_ci.csv` has one row per internal node: `node_id`, `label`,
  `time`, `se`, `lower`, `upper`, `width`, `se_reliable`, `notes`.  Rows where
  the interval does **not** come from the msz236 eq (7) branch-variance
  decomposition carry `se_reliable = False` and an explanatory `notes` entry
  (a node pinned at a calibration bound, an unsatisfied constraint, a
  lower-bound clip at 0).  The count is mirrored as
  `n_nodes_se_not_reliable` in the report JSON.

### `branch_var` / `n_sites` / `seq_length` — sampling variance vS(b)
* Default: `None` for all three (`ci.py:confidence_interval`).  `None` is not
  a variance value, it means "this source was not supplied"; zero variance
  is what the fallback produces when nothing is supplied, not a default.
* Behaviour — the first supplied source wins:
  1. explicit `branch_var` dict `{child node id: variance}`;
  2. Poisson approximation `vS(b) = b/L` with `L = seq_length`, otherwise
     `L = n_sites` (`seq_length` is tested first, so passing both uses
     `seq_length`);
  3. if all three are `None` (the default call) `vS(b) = 0`, so the interval
     captures only rate-heterogeneity variance (the msz236 simulation
     protocol).
* CLI: `--branch-var TSV` (`node_id`, `var` columns) and `--n-sites L`;
  `seq_length` has no CLI flag (API only).
* The source actually used is recorded in `<prefix>_report.json` as
  `confidence_interval.parameters.v_s_source`.  The Poisson approximation is
  documented as approximate (msz236 §sampling variance).

### `level`
* Default `0.95`.  CI bounds are truncated at hard calibration bounds, and a
  symmetric interval whose lower end would fall below 0 is clipped at 0 (a
  divergence time cannot be negative; `notes` records the clip).

### Rate-heterogeneity component `RV(R)` (not user-settable)
* `RV(R) = max(Vobs(R) - SV(R)/N, 0)` on the **per-lineage** scale of msz236
  eq (9): `SV(R)` is a sum over the `N` lineages, so it is divided by `N`
  before being subtracted from the per-lineage average `Vobs(R)`.  Comparing
  the raw sum against the average clamped `RV(R)` to 0 on every real tree and
  silently dropped the whole heterogeneity component.
* `Vobs(R)`, `SV(R)` and `RV(R)` are recorded in the report JSON as
  `Vobs_R`, `SV_R` and `RV_R`; the method string additionally echoes the
  sampling-variance source, the reproduced `method` and `mean`.
* When `SV(R)/N > Vobs(R)` (a sampling-dominated tree: too few sites for the
  observed rate spread) `RV(R)` clamps to 0, `v(b) = vS(b)`, and a warning
  says so — the interval then reflects sampling noise alone.
* A node whose calibrated age is a fixed point of its own constraint (on a
  hard bound, or the midpoint of the feasible interval) has an exactly zero
  derivative through eq (7).  Rather than reporting "zero uncertainty", such
  a node receives the uncertainty its calibration asserts (effective bounds,
  density quantiles, or the user `min`/`max` pair) and is flagged
  `se_reliable = False`.

## CorrTest (`corrtest`)

### `sister_resample`
* Default `0`.  For trees with fewer than ~50 tips use >= 50 resamples of
  the sister-pair assignment; the mean Spearman rho across resamples is
  used (R3F behaviour).

### `anchor_node`, `anchor_time`
* Defaults: `anchor_node=None`, `anchor_time=0.0` (`corrtest.py:corrtest`).
* `anchor_node=None` (the default) skips the absolute-rate step entirely;
  `anchor_time` is then unused.  Node ids are 1-based, so `0` is *not* the
  "no anchor" sentinel here (unlike R3F, where `anchor.node = 0` means "no
  anchor"): `anchor_node=0` raises `ValueError`, and an id that is not an
  internal node of the ingroup tree raises `KeyError`.
* CLI: `--anchor-node INT` (unset = `None`) and `--anchor-time FLOAT`
  (default `0.0`).
* Convert relative rates to absolute rates around a known node age; the
  mean/SD of absolute rates are appended to the report.

### Reference value and precision (not user-settable)
* Gate G6 compares against R3F's CorrScore on the bundled benchmark.  R3F
  prints only `format(score, digits = 5)`, so
  `data/golden/r3f/example_sr0_corrtest.txt` reads `0.9996`; the unrounded
  reference is `0.99959505940763882`
  (`data/golden/r3f/example_corrtest_full_precision.txt`, captured by
  `reproduce/golden_r3f_fullprec.R`).  G6 asserts bit-for-bit equality with
  the unrounded value, plus the identity between every rate CorrTest consumes
  and the engine's own per-node rate.
* With `--sister-resample > 0` the comparison is at the P-band level only:
  R3F's resampling draws come from R's RNG, which OpenRelTime cannot
  reproduce bit-for-bit even under the same seed.

## ddBD (`ddbd`)

### `sampling_frac`
* `None` (default) reports the rho the fit found.  A value is the sampling
  fraction to **report**: as in R3F (`ddbd.R` lines 491-505) the likelihood
  optimisation always leaves rho free with bounds `(0, 1)`, and the supplied
  value merely replaces the reported fraction — birth and death rates remain
  the optimum of the free-rho fit.  So `sampling_frac` does **not** restrict
  the fit, it overrides one reported number; the free-rho optimum stays
  visible as `sampling_frac_fitted` in the report parameters.
* Valid range `(0, 1]`; anything else raises `ValueError`.  Note that `0`,
  R3F's sentinel for "estimate it", is not accepted — use `None`.

### `anchor_node`, `anchor_time`
* Defaults: `anchor_node=None`, `anchor_time=1.0` (`ddbd.py:ddbd`).  Note the
  asymmetry with `corrtest`, where `anchor_time` defaults to `0.0`.
* CLI: `--anchor-node INT` (unset = `None`), `--anchor-time FLOAT`
  (default `1.0`).
* `anchor_node=None` (the default) leaves the relative times unscaled (scale
  factor 1), so the reported rates are per relative-time unit.  As in
  `corrtest`, `anchor_node=0` raises `ValueError` (node ids are 1-based; `0`
  is R3F's "no anchor" sentinel, not a valid id here) and an id that is not
  an internal node raises `KeyError`.
* Scale relative times so the anchor node has the given age.  Keep the
  maximum anchor time <= 10 for numerical stability (R3F guidance).

### `measure`
* `"SSE"` (default) | `"KL"`: initial-value grid selection score.  Grid:
  birth 1.1-10.1 step 1, death 1-10 step 1, sampling fraction
  {0.001, 0.01, 0.1, 0.5, 0.9}, filtered to birth >= death; optimisation by
  L-BFGS-B on the birth-death likelihood, trying starts in score order (max
  50).

## BLB (`blb`)

### `gamma`
* Default `0.7`, allowed range [0.6, 0.9]: subsample size `ceil(L**gamma)`.
* The default sits *inside* the admissible range of the RelTime-JA sampler
  (`lb_sampler.R` line 4); it is **not** the value the published RelTime-JA run
  used.  That run used `g = 0.78` with `20` subsamples x `20` replicates
  (`lbs_codelines.r` line 14), whereas OpenRelTime's defaults are `0.7` and
  `10 x 10`.  The defaults therefore reproduce the *method*, not that specific
  published protocol (see `docs/validation.md`).

### `n_subsample`, `n_replicate`
* Defaults 10 and 10 (100 IQ-TREE runs).  `seed` controls all resampling,
  including the per-replicate IQ-TREE `-seed`.  These are convenience defaults,
  not the RelTime-JA `20 x 20` protocol; raise them (with `gamma = 0.78`) to
  match a published BLB budget.

### `iqtree_args`
* Type: list of strings, e.g. `["-m", "GTR+G", "-B", "1000"]`.
* Default: `None` (`bootstrap.py:blb`), treated as "no extra arguments"
  (`iqtree_args or []`); an empty list is equivalent.
* CLI: `--iqtree-arg` is repeatable (`--iqtree-arg=-m --iqtree-arg=GTR+G`)
  and unset yields an empty list.
* Forwarded verbatim to IQ-TREE.  Tokens are validated against a strict
  allow-list and the process is started without a shell.

### `workdir` / `output_dir` (BLB)
* CLI `--workdir PATH` and `-o/--output-dir PATH` (default `ORT_blb`).
  Intermediate replicate alignments and trees are written under `workdir`
  (default: a temporary directory that is removed afterwards);
  `blb_summary.csv` and `blb_replicates.csv` are written into `output_dir`,
  which is created if missing.  Keep `--workdir` when an analysis must be
  audited replicate by replicate.

## MEGA-CC bridge (`megacc`)

### `megacc`
* CLI `--megacc PATH`, required: the path to the `megacc` executable.  It is
  not installed by CI and never downloaded — see `docs/validation.md`.

### `calibrations`
* CLI `-c/--calibrations PATH`, optional: the same `calibrations.tsv` format
  as `calibrate`.  Without it the bridge runs RelTime from branch lengths
  alone (`ppRelTimeBLens = true`).

### `mao`
* CLI `--mao PATH`, optional: use a hand-written `.mao` analysis-description
  file instead of the generated minimal one.  The generated file is written
  into `workdir` and printed by the run.

### `workdir` / output
* CLI `--workdir PATH`: where the `.mao` and megacc's stdout/stderr land.
* Results are written to `<prefix>_megacc_times.csv`, an R3F-compatible node
  time table whose `source` column records where every row came from.  Every
  argument is checked against an allow-list and megacc is started without a
  shell; if megacc produces nothing parseable the run raises rather than
  silently substituting OpenRelTime's own estimates.

## Numerical constants (advanced)

These are not exposed as options; they are listed because they determine
reproducible behaviour and mirror R3F.  All live in
`openreltime/_constants.py` with provenance comments.

| Constant | Value | Meaning | Source |
|---|---|---|---|
| `EPS_R3F` | 1e-19 | zero-branch substitute used exactly where R3F writes `10e-20` | `rrf_times.R` |
| `EPS_R3F_SMALL` | 1e-20 | zero-branch substitute used where R3F writes `1e-20` | `rrf_times.R` |
| `ZERO_BRLEN_WARN_FRACTION` | 0.10 | share of zero-length branches that triggers a warning | `rrf_times.R:67` |
| `RATE_RATIO_THRESHOLD` | 20.0 | rate-ratio guard threshold used by `rrf_times` and `rates-times`; the single default registered for `rate_ratio_threshold` | `rrf_times.R:354-356` |
| `RATE_RATIO_GUARD_DISABLED` | `None` | the sentinel that switches the guard off; the guard does not run in a rates-only pass | `rrf_rates.R:250-336` vs `rrf_times.R:353-373` |
| `CALIBRATION_MAX_ITER` | 100 | rounds of violation propagation before a conflicting calibration set is reported | msz236 / R3F |
| `EFFECTIVE_BOUNDS_REPLICATES` | 10 000 | replicates of the effective-bounds procedure | msz236 |
