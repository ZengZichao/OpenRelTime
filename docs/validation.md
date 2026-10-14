# Validation gates, coverage gaps and tested platforms

**Language: English** · [中文](validation-zh.md)

This file is the single place that states what OpenRelTime's automated suite
actually proves, what it does not prove, and which platforms the CI exercises.
Nothing here is softened, and every
"declared rather than demonstrated" claim is called out as such.

Command to reproduce the engine suite (needs the `dev` extra):

```bash
python -m pytest            # collection follows [tool.pytest.ini_options]: tests/
```

The desktop GUI is out of scope here: it ships as the independent
`OpenRelTime-Studio` project and is validated there, not by this suite.

## 1. Gate map (G1-G9) — what each gate asserts

`where` names the test file that carries the assertion.  "Asserted" means a
`pytest` assertion exists; "measured" means a number is computed/printed but no
pass/fail bound is enforced on it.

| Gate | Subject | Asserted condition | Where |
|---|---|---|---|
| G1 | Hand-computed RRF (msy044) | 4-taxon geometric + arithmetic node values match hand-worked constants to 1e-12; node-id/ape-cladewise contract | `tests/test_g1_hand_computed.py` |
| G2a | R3F rate regression (274-tip tree) | slope >= 0.999 **and** median relative deviation <= 1e-6 | `tests/test_regression_golden.py::test_g2a_rate_regression` |
| G2b | R3F time regression (274-tip tree) | slope >= 0.999 **and** all clade times matched (clade-keyed) | `tests/test_regression_golden.py::test_g2b_time_regression` |
| G2c | NEXUS `[&rate=...]` round-trip | read *and* write directions preserve rate annotations (the R3F golden NEXUS keeps its rates) | `tests/test_calibration_ci_edges.py::test_nexus_roundtrip_with_rate_annotations`, `tests/test_regressions_treeio_peripherals.py::test_r3f_annotated_golden_nexus_keeps_its_rates` |
| G3 | CLI/API equivalence | API and CLI share one code path (`report.write()`), so every artefact a subcommand writes — node tables, annotated timetrees, the calibrated and effective-bound CSVs — is **byte-identical** between the two entry points, and the command line is byte-identical with itself across two runs at a fixed seed. The JSON report is compared as a document, not as bytes: the CLI records tree-loading provenance (`resolve_polytomy`, `seed`, resolved input paths) that the API has no reason to invent. CLI options reach the engine — including `--seed`, which must reach the polytomy resolver and the report; CLI errors are clean, not tracebacks. **All byte-identity here is same-machine, same-seed** — it is *not* a cross-platform claim | `e2e_tests/cases/test_05_cli_api_parity.py` (`test_cli_and_api_artefacts_are_byte_identical_on_one_machine`, `test_api_write_matches_cli_rates_times`), `tests/test_regressions_engine.py` (`test_*_cli_*`, `test_report_json_merges_across_analyses`), `tests/test_regressions_treeio_peripherals.py` (`test_cli_*`, incl. `test_cli_seed_reaches_read_tree`, `test_cli_random_resolution_is_seed_reproducible` and `test_cli_ci_restores_the_recorded_polytomy_seed`, `test_branch_var_help_matches_the_header_the_parser_requires`, `test_rates_help_and_report_admit_the_guard_is_ignored`) |
| G4 | Input validation | negative branch length, duplicate tip names, polytomy (reject / `random` resolve, on the API and through the CLI), <3 ingroup tips, and "no internal node survives as a tip" guards raise | `tests/test_calibration_ci_edges.py` (`test_negative_*`, `test_duplicate_*`, `test_polytomy_*`, `test_two_tip_ingroup_rejected`), `tests/test_regressions_treeio_peripherals.py::test_validation_rejects_an_internals_degenerated_into_a_tip`, `tests/test_regressions_treeio_peripherals.py::test_cli_resolve_error_still_refuses_the_polytomy` |
| G5 | Outgroup rooting / taxonomy | rooting + pruning reproduces `ape::root(resolve.root=TRUE)` + `drop.tip` in topology and root-to-tip depths (compared as canonicalised `(tip, depth)` pairs, not as raw strings); complete-clade check (error/warn); group-name outgroups and calibrations; taxon-table parsing | `tests/test_regression_golden.py::test_rooted_topology_matches_ape`, `tests/test_taxonomy_monophyly.py` |
| G6 | CorrTest | CorrScore agrees with R3F's **unrounded** reference 0.99959505940763882 (the stored golden `example_sr0_corrtest.txt` only keeps 5 digits, 0.9996); lag-2/lag-3 decay features are non-degenerate; anchor scaling; every rate CorrTest consumes equals the engine's own per-node rate for **all 541 nodes of the benchmark** (270 internal + 271 tip rates, within 1e-9 absolute; on this tree the agreement is exact) | `tests/test_regression_golden.py::test_corrttest_matches_r3f_score`, `::test_corrttest_lag_decay_features_are_not_degenerate`, `tests/test_regressions_engine.py` (`test_corrtest_anchor_*`) |
| G7 | Calibration + analytic CI | bounds hit exactly; joint constraints; conflicting constraints raise — including a converged solution that is **not nested in the tree** (a negative age or a descendant older than its ancestor is rejected by name, not warned about); one-sided bound sets take their finite edge and are **invariant to the time unit** of the supplied ages; effective-bounds + taxon-set calibration, with the density row's effective bounds contained in **its own** density quantiles; CI identity at zero variance, widening with sampling variance, hard-bound truncation; a node sitting on its bound (to within the pinning tolerance) is **not** labelled "constraint not satisfied"; branch lengths outside the substitutions-per-site range are warned about once without changing the estimator; **bound-less calibration rows rejected; `RV(R)` positive and saturating; pinned nodes get calibration uncertainty instead of zero width; CI derivatives follow the calibrated `mean`/`method`** | `tests/test_calibration_ci_edges.py` (`test_bounds_*`, `test_two_calibrations_*`, `test_conflicting_*`, `test_effective_method_rejects_the_non_nested_format_demo`, `test_time_order_check_names_negative_and_inverted_nodes`, `test_one_sided_*`, `test_calibration_by_taxon_set_and_effective`, `test_node_on_its_bound_is_not_reported_as_unsatisfied`, `test_off_scale_branch_lengths_are_warned_once`, `test_ci_*`), `tests/test_regressions_calibration_ci.py` (`test_calibration_without_any_bound_is_rejected`, `test_rate_heterogeneity_*`, `test_pinned_nodes_*`, `test_arithmetic_point_estimate_*`, `test_rrf_and_implied_rates_*`, `test_clamp_warning_*`, `test_divergent_trial_ids_are_detected`) |
| G8 | ddBD | fitted (birth, death, sampling) match the R3F run within optimiser tolerance; sampling-fraction semantics; grid order; L-BFGS-B abnormal-termination handling | `tests/test_regression_golden.py::test_ddbd_close_to_r3f`, `tests/test_regressions_engine.py::test_ddbd_*`, `tests/test_regressions_treeio_peripherals.py::test_ddbd_*` |
| G9 | Runtime / scalability | **measured, not asserted.**  Per-node timing and the 1,000-tip < 1 s figure come from a single Apple M-series laptop; no CI job asserts any runtime bound.  The largest tree the suite builds is a synthetic **1,200-tip pectinate** chronogram (`tests/test_regressions_engine.py::test_deep_tree_full_pipeline`), which is run only to prove the parser, traverser and serialiser survive deep recursion — it is never timed — and no CI job runs a 5,000+-taxon tree at all | — (no automated gate) |

**Notes on gate strength (not hidden)**

* **G2a (rates)**: the *gate* asserts median relative deviation <= 1e-6 and
  never checks the maximum, so the rate maximum below is an observation rather
  than an assertion.  **G2b (times)**: the *gate* asserts the per-node *maximum*
  relative deviation <= 1e-6 (over the 541 shared clades, alongside slope >=
  0.999), so it does bound the worst node.  Measured on the 274-tip benchmark,
  the medians are ~6.8e-16 (rates) and ~8.7e-16 (times) and the maxima
  ~4.8e-15 / 3.6e-14; the rate maximum carries no assertion, whereas bounding
  the time maximum is exactly what G2b does.
* **G1**: the arithmetic-mean convention is only exercised on shallow
  (2-layer) trees, where geometric and arithmetic folding coincide; deeper
  arithmetic behaviour is not hand-gated here.
* **G3**: "byte-identical" is a same-machine, same-seed property (CSV locale,
  line endings, JSON key order and zip metadata can differ across platforms);
  it is deliberately not stated as a cross-platform guarantee.  The CLI/API
  byte test compares the *data* artefacts as bytes and the JSON report as a
  document, because the CLI records tree-loading provenance the API does not.
* **G6/G8**: these are regression checks against R3F's stored output on the one
  bundled benchmark tree; they do not establish agreement on other topologies.

## 2. Shipped modules and their validation status

`openreltime/` ships these Python modules.  Status is what the suite actually
proves, not what the module advertises.

| Module | Status | Basis |
|---|---|---|
| `rrf`, `times`, `rates` | **validated** vs R3F + hand-computed | G1, G2a, G2b |
| `tree`, `treeio` | **validated** (topology, rooting, node-id contract, NEXUS) | G2c, G4, G5 |
| `calibrate` | **validated** (bounds, conflicts, effective) | G7 |
| `ci` | **validated** (structure: identity/widening/truncation) | G7 |
| `corrtest` | **validated** vs R3F score on the benchmark | G6 |
| `ddbd` | **validated** vs R3F params on the benchmark | G8 |
| `report` | **validated** (write flags, JSON merge) | G3 |
| `taxonomy` | **validated** (monophyly, group resolution, table parsing) | G5 (`test_taxonomy_monophyly.py`) |
| `table` | **validated** vs R3F `tree2table()` golden | `test_tree2table_matches_golden` |
| `_external`, `_constants`, `cli` | infrastructure; behaviour smoke-tested via CLI/G4 tests, no standalone reference gate | — |
| **`bootstrap` (BLB)** | **shipped but NOT validated numerically.**  CI never installs IQ-TREE; the one test stubs the executable to check argv/`-seed` plumbing and degenerate-replicate guards.  The resampling/aggregate *distribution* is not compared to any reference | no reference gate |
| **`megacc`** | **shipped but NOT validated.**  The bridge needs a real `megacc` binary, which is not installed or run by CI; parse tests feed *synthetic* outputs.  MEGA-CC agreement is established against the published R3F/MEGA results and the shipped golden files, not by a head-to-head run | no reference gate |
| **`viz`** | **shipped but NOT validated against a baseline.**  Tests assert structural properties (log colour ramp; ddBD density integrates to 1 on the plotted axis); there is no golden-image comparison | no reference gate |

`bootstrap`, `megacc` and `viz` are therefore described elsewhere in this
manual as **shipped but not validated here**; `taxonomy` and `table` are covered
by the gates listed above.

## 3. Platforms and Python versions — tested vs claimed

From `.github/workflows/ci.yml`:

| Axis | Tested by CI | Merely claimed elsewhere |
|---|---|---|
| OS | `ubuntu-latest`, `macos-latest` (matrix `test` job) plus a separate `windows-latest` job (`test-windows`) running both the engine and the end-to-end suites | Windows support was **declared but undemonstrated** up to 0.1.0; the job is declared for 0.1.1 and its first run is the next push.  See §5 for the one Windows-specific code limitation the suites do not cover |
| Python | **Ubuntu:** 3.10, 3.11, 3.12, 3.13. **macOS: 3.13 only** — `ci.yml` excludes `macos-latest` × 3.10/3.11/3.12, so the five matrix entries are 4 Ubuntu + 1 macOS, not 4 × 2 | any "3.10–3.13 on Linux **and** macOS" wording is wrong for the engine (the separate Studio project's CI does run 4 × 2) |
| Desktop GUI | **not tested here** — the GUI is the separate `OpenRelTime-Studio` distribution | its headless/interactive coverage lives in that project's own suite and CI |
| External tools | **not installed** — `blb` (IQ-TREE) and `megacc` (MEGA-CC) paths are never executed in CI | — |
| Runtime scale | the golden/regression benchmark is the 274-tip tree; the largest tree the suite builds at all is a synthetic 1,200-tip pectinate chronogram, run only for parser/traverser/serialiser stack depth (`tests/test_regressions_engine.py::test_deep_tree_full_pipeline`) and never timed | the "1,000 tips < 1 s" figure is from one laptop, **not** a CI assertion (G9), and no timed benchmark above 274 tips runs in CI |

The matrix above is what `.github/workflows/ci.yml` declares and executes on
every push and pull request; its run logs are published on the repository's
Actions page.  A platform outside that matrix is a declared target, not a
demonstrated one.

## 4. External tools — citation and assumed version

**IQ-TREE** (invoked by `openreltime/bootstrap.py` for BLB; not pinned anywhere
in the code or CI):

* Citation: Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Wood MR, von
  Haeseler A, Lanfear R. *IQ-TREE 2: New Models and Efficient Methods for
  Phylogenetic Inference in the Genomic Era.* Mol Biol Evol 37:1534-1538
  (2020), doi:10.1093/molbev/msaa156.
* Version actually assumed: **none is pinned or installed by CI.**  The flags
  used (`-s`, `-pre`, `-nt AUTO`, optional `-seed`, `-m`) and the
  `<prefix>.treefile` output are the IQ-TREE **2.x** command-line convention
  (retained by IQ-TREE 3.x).  Development/figures used the IQ-TREE 2.x/3.x
  line; treat any specific minor version as the operator's choice, and record
  `iqtree --version` output alongside a BLB run.

**MEGA-CC / MEGA equivalence boundary (one sentence)**: OpenRelTime implements
only the analytical RRF layer, whereas MEGA's RelTime applies an additional
likelihood-ratio-test rate-constraint refinement that OpenRelTime does not, so
the two are expected to be highly correlated (~0.98-0.99) but not equal; MEGA-CC
agreement is established against the published R3F/MEGA results and the shipped
golden files, not by a new head-to-head run.

## 5. Known behavioural limitations the gates do NOT catch

These are restated here so that the green gates are not over-read:

* `--outgroup` is optional in the CLI; an unrooted input is accepted and its
  root position silently drives every age.  R3F requires an outgroup.  No gate
  covers this: `openreltime times -i tree.nwk` without an outgroup runs to
  completion with no warning, so the caller must confirm the input is rooted on
  the intended clade.
* The absolute width of an analytical interval depends on the unit of the input
  branch lengths.  msz236 eq (14) spreads `RV(R)` over branches by the ratio
  `b_j^2 / sum(b^2)`, which is scale-free, but `vS(b) = b / L` and the
  derivatives `dt/db` are not, so the same tree expressed in substitutions per
  100 sites rather than per site is qualified differently.  The estimator is the
  published one and is deliberately not altered; `openreltime ci` now warns once
  when the longest branch falls outside 0.02-5 substitutions per site, but the
  number it prints still depends on what you fed it.
* At the v0.1.1 tag the two external-tool bridges still could not start on
  Windows: the command-token blacklist of `openreltime/_external.py`
  (`_FORBIDDEN`) rejected backslashes and double quotes, which are legitimate
  in a Windows path, so `resolve_executable()` refused the `which("iqtree")`
  result and `megacc`'s `_safe_path()` refused every absolute
  `.mao`/tree/output path.  The blacklist on `main` no longer contains the
  backslash — the fix landed after the v0.1.1 tag, in the first
  `test-windows` run, and ships with the next release — and with `shell=False`
  a backslash in a list-mode token is inert data, so Windows paths validate; a
  unit gate pins a `C:\Program Files\...`-style token as acceptable.
  End-to-end runs of `blb` and `megacc` remain *declared, not demonstrated* on
  Windows, because CI installs neither IQ-TREE nor MEGA-CC; the analytical
  core and the CLI are exercised there by the `test-windows` job.

## 6. Calibration and CI behaviours under regression gate

Each behaviour below carries an assertion in
`tests/test_regressions_calibration_ci.py` that fails if it regresses, so none
of them is a coverage gap; the user-visible description of the behaviour is what
`docs/parameters.md` and `docs/methods.md` state.

| Behaviour asserted | Test |
|---|---|
| a calibration row carrying neither bound is rejected when the table is parsed, naming the offending line; a non-finite feasible factor raises | `test_calibration_without_any_bound_is_rejected`, `test_parse_calibrations_names_the_offending_line`, `test_feasible_factor_raises_instead_of_returning_inf` |
| `RV(R) = max(Vobs(R) - SV(R)/N, 0)` is evaluated on the per-lineage scale, so the msz236 heterogeneity component is positive in the few-sites regime and saturates at `Vobs(R)` as sites grow; a warning says when the clamp still bites | `test_rate_heterogeneity_variance_is_positive_at_low_site_counts`, `test_rate_heterogeneity_saturates_at_Vobs_as_sites_grow`, `test_clamp_warning_*` |
| a node pinned at a hard bound receives the uncertainty its own calibration asserts rather than a zero-width interval, and is flagged `se_reliable = False` with a `notes` explanation | `test_pinned_nodes_get_the_calibration_uncertainty`, `test_midpoint_fixed_calibration_is_not_zero_width`, `test_unsatisfied_constraint_never_inverts_the_interval` |
| the differentiation pipeline reproduces the calibrated `mean`, the guard threshold and, for `method="effective"`, the effective bounds | `test_arithmetic_point_estimate_gets_arithmetic_derivatives`, `test_effective_result_is_differentiated_along_its_effective_bounds`, `test_divergent_trial_ids_are_detected` |
| post-calibration lineage rates and RRF rates are named apart, and the variance decomposition uses the RRF-consistent definition | `test_rrf_and_implied_rates_are_named_apart`, `test_variance_uses_the_rrf_rates_and_drops_lineages_without_duration`, `test_non_finite_rates_are_dropped_from_the_average` |
| the repair loop's per-round nested check is followed by a total re-check: a converged vector with a negative age or a descendant older than its ancestor raises and names the nodes, so the shipped format-illustration table is rejected under `method="effective"` instead of returned with warnings, while the feasible example still solves to `f = 371.59973111298007` with zero violations | `tests/test_calibration_ci_edges.py::test_effective_method_rejects_the_non_nested_format_demo`, `::test_time_order_check_names_negative_and_inverted_nodes`, `::test_feasible_example_calibrations_still_solve_cleanly`, `e2e_tests/cases/test_02_calibration_ci.py::test_calibrate_effective_rejects_the_format_demo` |
| a one-sided bound set takes its finite edge on the relative scale, so the same physical constraint gives ages 1000x apart when stated in Ma and in Ga, and a bound younger than the uncalibrated solution is no longer left inert | `tests/test_calibration_ci_edges.py::test_one_sided_min_bound_below_the_uncalibrated_age_is_not_inert`, `::test_one_sided_max_bound_above_the_uncalibrated_age_does_not_inflate`, `::test_one_sided_bounds_are_invariant_to_the_time_unit` |
| the effective bounds of a density row stay inside that density's own 2.5/97.5 quantiles, and a node sitting on its bound is not written into the `notes` column as a constraint the solver failed to reach | `tests/test_calibration_ci_edges.py::test_calibration_by_taxon_set_and_effective`, `::test_node_on_its_bound_is_not_reported_as_unsatisfied`, `::test_genuinely_unsatisfied_constraint_is_still_flagged` |
| branch lengths outside the substitutions-per-site range earn one warning without changing the estimator, because eq (14) makes the absolute half-width depend on the unit of the input | `tests/test_calibration_ci_edges.py::test_off_scale_branch_lengths_are_warned_once` |

The corresponding user-visible documentation lives in `docs/parameters.md`
(sections *Calibration*, *Confidence intervals*) and `docs/methods.md` (§4).

## 7. Precision of the CorrTest reference value

`data/golden/r3f/example_sr0_corrtest.txt` stores only five significant digits
because R3F prints `format(score, digits = 5)`. The unrounded reference double is
in `data/golden/r3f/example_corrtest_full_precision.txt` (see
`reproduce/golden_r3f_fullprec.R` for how it was captured without touching the
numerical path). Gate G6 asserts equality with that unrounded value, which is a
bit-for-bit identity, and separately asserts that every rate CorrTest consumes
matches the engine's own per-node rate for all 541 nodes of the benchmark — the
270 internal rates and the 271 tip rates, which agree exactly on this tree. The
golden directory itself was
re-generated in the recorded R environment and reproduced eleven of twelve files
byte for byte; the NEXUS timetree differs only in ape's embedded write-time stamp.
