# Changelog

**Language: English** · [中文](CHANGELOG-zh.md)

OpenRelTime follows Semantic Versioning.  This file records the contents of
each release, newest first.

## 0.1.1 — 2026-10-08

Bug-fix and hardening release prepared during pre-submission review of the
accompanying *Bioinformatics* Application Note.  Every item below was found by
auditing the manuscript's claims against the code and then reproducing the
defect.  Three of them could silently produce a wrong answer, so they are
treated as corrections rather than refinements.

### Fixed

- **Unsatisfiable nested calibrations now raise.**  A bound that required a
  descendant node to be older than its ancestor previously appended a warning
  and returned a solution that violated the requested `min_bound`, while also
  abandoning repair of every remaining calibration row.  It now raises
  `ValueError` naming the node and both bounds, matching the existing
  root-cannot-be-younger path.
- **That guarantee now covers the converged solution, not just one round.**
  The nested check above examined only the worst violator of the round being
  repaired, so a solution could still finish non-monotonic and be returned with
  warnings.  Measured on `data/examples/example_calibrations.tsv` under
  `method='effective'` (`n_effective=10000`, `seed=20260906`):
  `time_factor` 4432.82, root 1162.86 Ma, node 307 at **-3998.56 Ma**, four
  parent-child inversions, four `NaN` implied rates, no error.  Every emitted
  node age is now re-checked for non-negativity and for the tree's own
  parent-child ordering, and a failure raises naming the offending node ids.
  `example_calibrations_feasible.tsv` is unaffected (`f = 371.59973111298007`,
  zero violations).
- **One-sided calibration bound sets take their finite edge.**  The feasible
  factor was floored at `1.0`, a constant of the unitless RRF relative scale
  rather than of the caller's time unit, so a minimum bound that asked for a
  younger tree than the uncalibrated solution was left silently inert (bound
  0.1 → `f = 1.0`) while a maximum bound above it *inflated* every age in the
  tree (bound 0.3 → `f = 1.2862`), and the same physical constraint behaved
  differently when stated in Ma and in Ga.  The edge is now used directly, which
  is what the documentation said all along and makes the result invariant to the
  unit.  Two-sided sets — every shipped example and benchmark file — are
  untouched, because they never reached that branch.
- **The CI `notes` column no longer calls a satisfied bound unsatisfied.**
  A node the solver had placed on its hard bound was labelled "calibration
  constraint not satisfied by the solver" whenever it missed the comparison by
  one ULP, writing a false statement into a user-facing CSV column, while the
  neighbouring pinned test already used the pinning tolerance.  Both now use it.
- **Off-scale branch lengths are announced.**  msz236 eq (14) is scale-free in
  the `b_j^2 / sum(b^2)` ratio but not in the resulting half-widths, so the
  absolute interval width depends on the unit of the input branch lengths
  (measured: the same tree gives widths of 0.426197, 0.042620 and 0.0042620 as
  its branches are scaled by 1, 1/10 and 1/100).  The estimator is the cited one
  and is deliberately not altered; `openreltime ci` now warns once per call when
  the longest branch falls outside 0.02-5 substitutions per site.
- **`corrtest(anchor_node=...)` no longer emits `inf`/`NaN`.**  Its
  `anchor_time` defaults to `0.0`, so the anchor scaling divided by zero and
  wrote `{"mean.rate": Infinity, "sd.rate": NaN}` into the JSON report — which
  is not valid strict JSON and is rejected by conforming parsers.  A
  non-positive `anchor_time` with an anchor node is now rejected up front, and
  reports are serialised with `allow_nan=False`.
- **Unknown density parameter keys are rejected instead of ignored.**
  `density='normal'` with `density_params='mu=5;sigma=1'` was accepted and
  sampled from the default `N(0, 1)`, silently applying the wrong prior.  Each
  family now has a key whitelist (`uniform`: min,max; `exponential`:
  offset,mean; `normal`: mean,sd; `lognormal`: offset,meanlog,sdlog) and
  unknown or missing keys raise.
- **Degenerate CorrTest inputs no longer report significance.**  An undefined
  Spearman correlation was mapped to `0.0` and fed to the fixed logistic model,
  so an all-zero-branch-length tree returned `CorrScore ≈ 1` with band
  `P < 0.001` — the strongest possible claim from zero information.  Such trees
  are now rejected.
- **`n_effective=1` is rejected** rather than silently running two replicates
  while reporting one.
- **`openreltime.<module>` name shadowing.**  The package re-exported functions
  whose names equalled their modules, so `openreltime.calibrate.Calibration`
  raised `AttributeError`.  Those attributes are now callable module objects:
  both `openreltime.calibrate(...)` and `openreltime.calibrate.Calibration`
  work.
- **`--seed` now reaches random polytomy resolution.**  `read_tree()` accepted
  a `seed`, but the CLI's tree-loading helper never passed one, so
  `--resolve random` was irreproducible and unrecorded on every subcommand
  except `blb`.  Every command that takes `--resolve` now also takes `--seed`
  and forwards it; the seed is recorded in the report (`monophyly`, which
  writes only a CSV and has no report block, is the documented exception), and
  `ci` re-reads the tree with the seed recorded by the calibration it qualifies.
- Division-by-zero guards on the normalisation maxima in `ddbd` and `corrtest`.

### Added

- `blb` and `megacc` now write JSON reports and record the version banner of
  the external tool they invoked, so a run can be attributed to a specific
  IQ-TREE or MEGA-CC build.  No release is pinned.
- A FigTree-readability test for the emitted NEXUS (rooted-tree token,
  annotation-free `TAXA` block, one `[&rate=...]` comment per node, numeric
  branch lengths, independent re-parse), plus an explicit test that plain
  `.nwk` output does **not** carry rate annotations.
- `data/examples/example_calibrations_feasible.tsv`: a runnable, mutually
  consistent example.  The pre-existing `example_calibrations.tsv` is a syntax
  illustration only — its density-only row makes it raise under the default
  `method='bounds'`, and its inconsistent bounds drive four lineages to
  non-positive durations.
- `tools/check_gates.py` cross-validates the gate table in
  `docs/validation.md` against collected test ids and `@pytest.mark.gate`
  marks; `tools/count_tests.py` reports per-suite counts.
- `@pytest.mark.gate(...)` marks across the unit suite, registered in
  `pyproject.toml`.
- `figures/`: the manuscript's four figure scripts and the two `-data.json`
  files they draw from, plus `figures/README.md`.  They are byte-identical to
  the sources that rendered the submitted artwork — re-running each script with
  `--outdir` reproduces all four published PNGs byte for byte — and
  `MANIFEST.in` ships them in the source distribution, so the Data-availability
  statement that names `figures/` is true of the deposit rather than of the
  authors' working copy.  Rendered artwork stays out of the distribution.

### Changed

- `docs/validation.md`, `docs/validation-zh.md` and `README.md` claimed CI runs
  Python 3.10–3.13 on both Ubuntu and macOS.  `ci.yml` excludes macOS for
  3.10/3.11/3.12, so the documents now state the real matrix: Linux
  3.10–3.13, macOS 3.13 only.
- `ZERO_BRLEN_WARN_FRACTION` and `MIN_INGROUP_TIPS` are now actually referenced
  by the code that uses them, instead of being duplicated as literals.
- The `ci.py` module docstring no longer claims intervals use calibrated
  absolute rates; they use the scale-invariant relative RRF rates, as the code
  has always done.
- Gate ids are machine-checkable and `tests/test_taxonomy_monophyly.py` is
  labelled with the gate the documentation assigns it.
- Documentation now matches the code and the tests where it did not: the G5
  row no longer claims byte-identity for ape-equivalent rooting (the assertion
  compares canonicalised `(tip, depth)` pairs), the G2a/G2b note is split
  because G2b does bound the per-node maximum, the FigTree NEXUS contract cites
  the real test id, the `seed` section names the stdlib `random.Random` used for
  polytomy binarisation instead of promising "always
  `numpy.random.Generator`", and `data/golden/r3f/example_RRF_rates.nwk` finally
  has the provenance row that file's own header promises for every data file.

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
