# Changelog

**Language: English** · [中文](CHANGELOG-zh.md)

OpenRelTime follows Semantic Versioning.  This file records the contents of
each release, newest first.

## 0.1.0 — 2026-09-29

First release of the OpenRelTime engine: a GPL-3 Python implementation of the
**relative rate framework (RRF)** of Tamura, Tao & Kumar (2018, *MBE*
35:1770-1782), the analytical core of RelTime in MEGA.  Relative lineage rates
and node times are estimated from a branch-length tree, calibrations are
converted to absolute times with analytical confidence intervals, and CorrTest,
ddBD and a bag-of-little-bootstraps pipeline are provided alongside.  Every
analysis is reachable both from Python and from the command line.

Requires Python >= 3.10 with numpy, scipy, pandas and click; no compiled
extensions.

### Analysis engine

- **RRF rate and time estimation** (`openreltime.rrf`, `rates`, `times`;
  `rrf_rates`, `rrf_times`, `rrf_rates_times`): closed-form solution of the
  local 3- and 4-lineage problems in postorder with a preorder adjustment pass,
  subtree folding into composite branch lengths, and the ingroup root anchored
  at average rate 1.
- **Two averaging conventions**: the geometric mean (msy044 eqs 28-42, the
  convention R3F and RelTime use; default) and the arithmetic mean (msy044
  eqs 1-27, the 2012 semantics), with the folding operator following the chosen
  convention.
- **Extreme-rate guard**: a node whose adjusted rate leaves the threshold band
  takes the time of the nearest non-exceeding ancestor (R3F semantics); default
  threshold `20.0`, disable with `--no-guard`, sensitivity count recorded as
  `n_rate_guarded_nodes`.
- **Numerical guards shared with R3F**: epsilon substitution at exactly the
  positions R3F writes it, and a warning once zero-length branches exceed 10%
  of the tree (all defaults registered in `openreltime/_constants.py`).

### Calibration and confidence intervals

- **Calibration to absolute times** (`openreltime.calibrate`): the `bounds`
  method (global time factor *f* as the midpoint of the feasible-interval
  intersection, proportional rate rescaling of the offending lineage propagated
  to descendants, capped at 100 iterations, conflicting sets reported with
  per-node diagnostics) and the `effective` method (msz236 density resampling,
  10,000 replicates by default).
- **Densities**: `uniform`, `exponential`, `normal` and `lognormal`, parsed
  strictly as `key=value` pairs (no expression evaluation), alongside hard
  `min_bound`/`max_bound` values in a tab-separated `calibrations.tsv`.
- **Node targeting** by internal node id or by the MRCA of a taxon set, whose
  entries may be tip labels or a group name.
- **Cancellable calibration**: `calibrate(..., progress=...)` returns
  `CalibrationCancelled` when the callback returns `False`.
- **Analytical confidence intervals** (`openreltime.ci`): the delta method of
  msz236 (eqs 7-14), with branch variance split into sampling and
  rate-heterogeneity components, partial derivatives obtained by central finite
  differences of the whole estimation pipeline, hard-bound truncation and
  clipping at 0.  Sampling variance comes from an explicit per-branch `var`
  table, a Poisson approximation `vS(b) = b/L`, or nothing at all.
- **One estimator throughout**: the differentiation reproduces the averaging
  convention, guard threshold and calibration bounds (including the effective
  bounds of an `effective` run) that produced the point estimate.
- **Reported reliability**: a calibrated node whose age is a fixed point of its
  own constraint is given the uncertainty its calibration asserts and flagged
  `se_reliable = False` with an explanatory `notes` entry; the effective-bounds
  method additionally reports empirical 2.5/97.5 node-age quantiles.

### Rate testing, tree priors and resampling

- **CorrTest** (`openreltime.corrtest`): the rate-autocorrelation statistic of
  Tao et al. (2019) with the published fixed-coefficient logistic model, P-band
  reporting, optional sister-pair resampling for small trees and optional
  absolute-rate anchoring.
- **ddBD** (`openreltime.ddbd`): birth-death speciation prior estimation (birth
  rate, death rate, sampling fraction) with an SSE/KL-scored initial-value grid
  and L-BFGS-B fitting, optional anchor-node/anchor-time scaling, and a
  plottable fitted density.
- **Bag of little bootstraps** (`openreltime.bootstrap`): two-level site
  resampling from an alignment, replicate trees inferred by an external IQ-TREE
  executable, per-clade aggregation of node ages into medians and 2.5/97.5
  percentiles with the fraction of replicates containing each clade, and the
  run seed forwarded to every IQ-TREE call.

### Tree input/output and taxonomy

- **In-house Newick and NEXUS I/O** (`openreltime.treeio`, `openreltime.tree`):
  `[&rate=...]` node annotations preserved on read and regenerated on write
  (R3F/FigTree dialect), comments and quoted labels handled without corrupting
  names.
- **Outgroup rooting and pruning** (`root_outgroup`) that reproduces the
  observed behaviour of `ape::root(resolve.root = TRUE)` + `ape::drop.tip`,
  including contraction of unary nodes and discarding the root edge, with a
  complete-clade check for multi-tip outgroups and a minimum of three ingroup
  tips.
- **ape-cladewise node numbering** (`PhyloNode.assign_ids`: tips 1..n, internal
  nodes in preorder, root = n+1), so node tables and `--anchor-node` values are
  directly comparable with R3F, ape and MEGA output.
- **Polytomy handling**: rejected by default, optionally resolved at random
  with zero-length branches.
- **Monophyly and classification tools** (`openreltime.taxonomy`):
  `check_monophyly`, `parse_taxon_table` (TSV/CSV tip-to-group tables) and
  `detect_groups_from_labels` (automatic `Genus_species` genus recognition);
  group names may also be used as outgroups and as calibration `taxon_set`
  values.
- **Node tables** (`openreltime.table.tree2table`): branch-length or node-age
  tables in R3F's layout.

### Outputs and reproducibility

- **Result objects** (`openreltime.report`) with `to_pandas()`, `to_json()` and
  `write(prefix)`; CSV, Newick, NEXUS and JSON outputs come from the same code
  path as the CLI, and `--r3f-compat` selects R3F's file and column names.
- **Reproducibility reports**: every analysis stage whose parameters can move
  the numbers (`rates`, `times`, `rates-times`, `calibrate`, `ci`, `corrtest`,
  `ddbd`) merges its parameters, warnings and random seed into
  `<prefix>_report.json`, one top-level key per run sharing that prefix; `ci`
  replays the settings its `calibrate` run recorded.
- **All randomness** (density sampling, BLB, sister resampling) routed through
  one seeded `numpy.random.Generator`.

### Interfaces

- **Python API** `import openreltime as ort` exposing the functions above plus
  `read_tree`, `parse_calibrations` and the result types.
- **Command line** `openreltime` (alias `ort`) with 11 subcommands mirroring
  the API: `rates`, `times`, `rates-times`, `calibrate`, `ci`, `corrtest`,
  `ddbd`, `tree2table`, `monophyly`, `blb` and `megacc`; readable error messages
  rather than tracebacks.
- **MEGA-CC bridge** (`openreltime.megacc`): `.mao` analysis-description
  generation, an external `megacc` runner, and a parser that returns an
  R3F-compatible node-time table in which a `source` column records the
  provenance of every row.
- **Optional plotting** (`openreltime.viz`, `[plot]` extra): rate-coloured
  timetrees on a logarithmic colour ramp, node-age interval plots, and ddBD
  densities normalised on the plotted axis.
- **External processes are started without a shell** and every forwarded
  argument is checked against a strict allow-list (`openreltime._external`).

### Validation and data

- **Reference-tool golden standard**: the R3F outputs on the 274-tip mammalian
  benchmark (`data/golden/r3f/`), the R script that produced them
  (`reproduce/golden_r3f.R`), the example inputs (`data/examples/`) and a
  per-file provenance record (`data/PROVENANCE.md`).
- **Automated gates** (`tests/`): G1 hand-computed 3-/4-taxon values for both
  conventions; G2a/G2b the R3F rate and time regressions on the 274-tip tree;
  G2c the NEXUS `[&rate=...]` round trip; G3 CLI/API equivalence; G4 input
  validation; G5 ape-equivalent rooting and taxonomy; G6 CorrTest against the
  unrounded R3F score; G7 calibration and CI behaviour; G8 ddBD parameters; G9
  runtime measurements (recorded, not asserted).  Numerical regression suites
  live in `tests/test_regressions_engine.py`,
  `tests/test_regressions_calibration_ci.py`,
  `tests/test_regressions_rrf_engine.py` and
  `tests/test_regressions_treeio_peripherals.py`.
- **End-to-end suite** (`e2e_tests/`): full CLI-and-API runs on the bundled
  example data against the same goldens, shipped with the source package.

### Documentation

- Bilingual manual (English and Chinese) throughout: `README.md`,
  `docs/usage-en.md`, `docs/parameters.md`, `docs/methods.md`,
  `docs/validation.md`, three tutorials in `docs/tutorials/`, the
  architecture-decision-record index in `docs/adr/`, `CHANGELOG.md`,
  `THIRD-PARTY-NOTICES.md` and `data/PROVENANCE.md`.

### Scope of this release

- **Analytical RRF only.**  MEGA's RelTime adds a likelihood-ratio-test
  rate-constraint refinement that this engine does not implement, so OpenRelTime
  agrees with MEGA in correlation (~0.98-0.99 on published benchmarks), not
  node-for-node; report correlations when comparing the two.
- **Reference agreement is scoped to the bundled benchmark.**  The CorrTest and
  ddBD gates are regression checks against R3F's stored output on the 274-tip
  tree and do not establish agreement on other topologies; the gate on the rate
  and time regressions asserts the slope and the *median* relative deviation,
  while maxima are measured and reported but not asserted.
- **`bootstrap`, `megacc` and `viz` ship without a numerical reference gate.**
  The BLB and MEGA-CC paths need external executables that the test suite does
  not install; their tests cover argument plumbing, output parsing and
  structural properties.  MEGA-CC agreement is established against the
  published R3F/MEGA results and the shipped golden files, not by a new
  head-to-head run.
- **Platforms.**  CI runs the suite on Linux and macOS across Python 3.10-3.13;
  Windows is a declared target that CI does not exercise.  Reported timings come
  from a single machine and the suite stops at the 274-tip golden tree.
- **The outgroup is optional.**  `--outgroup` may be omitted; the tree is then
  used exactly as read and its root position determines every age without a
  warning, so confirm the input is rooted on the intended clade (R3F requires an
  outgroup).
