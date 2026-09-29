# OpenRelTime End-to-End Test Plan

**Language: English** · [中文](TEST_PLAN-zh.md)

## 1. Objectives

This end-to-end (E2E) test plan aims to automate validation of the complete OpenRelTime analysis workflow using **real example data**. It covers the full spectrum from the command-line interface (CLI) to the Python API, from core algorithms to external-tool bridges, and from input validation to output reproducibility. All test files live in the repository's `e2e_tests/` directory and are shipped with the source package so users can reproduce the verification in any environment.

## 2. Scope

| Functional Area | Coverage |
|-----------------|----------|
| RRF core | `rates` / `times` CLI, golden regression, output files (CSV, NEXUS, Newick), normalization, arithmetic mean, guard options, NEXUS round-trip |
| Calibration & CI | `calibrate` (bounds/effective), `ci` (default, branch-var, n-sites), API vs CLI parity, conflicting calibrations, density-only constraints, pinned-node flagging |
| Correlation / birth-death | `corrtest` golden regression, full-precision API, resampling; `ddbd` golden, anchor, plot, API/CLI parity |
| Tree table & monophyly | `tree2table` branch lengths and node ages; `monophyly` auto/table/ad-hoc groups |
| CLI/API parity | API-written outputs match CLI, `_report.json` merging, parameter recording |
| Input validation | Negative branch lengths, duplicate tips, polytomies, incomplete outgroups, empty calibration rows, density-only vs method mismatch, missing input files |
| External tools | BLB (IQ-TREE); MEGA-CC (`megacc`) bridge |

## 3. Test Data

All data reside under `e2e_tests/data/` and never leave the working directory.

| File | Source / Generation | Purpose |
|------|---------------------|---------|
| `example.nwk` | `data/examples/example.nwk` | Main phylogeny |
| `example_taxa.tsv` | `data/examples/example_taxa.tsv` | Taxon grouping table |
| `example_calibrations.tsv` | `data/examples/example_calibrations.tsv` | Calibration constraints |
| `small.fasta` | Randomly generated (10 taxa × 300 bp) | BLB / IQ-TREE smoke test |
| `example.fasta` | Randomly generated (274 taxa) | Optional larger alignment sample |
| `golden/` | `data/golden/r3f/` | Regression baseline |

## 4. Environment

- **Operating system**: Linux x86_64
- **Python**: ≥3.10
- **Core dependencies**: numpy, scipy, pandas, click
- **Optional dependencies**: matplotlib, pytest
- **External tools**: IQ-TREE (BLB tests), MEGA-CC `megacc` (megacc bridge tests)
- **Installation**:
  ```bash
  pip install --editable ".[dev,plot]"
  micromamba install -p /path/to/openreltime/env -c conda-forge -c bioconda iqtree
  ```

## 5. Execution Strategy

- Use `pytest` as the test framework; test files are named `test_*.py`.
- Invoke `python -m openreltime.cli` via `subprocess.run` to exercise the real CLI entry point.
- External-tool tests use `skipif` to detect executables; they are skipped cleanly when absent.
- Each test uses an isolated `tmp_out` directory to avoid cross-test contamination.

## 6. Test Case Design and Acceptance Criteria

### 6.1 RRF core (`test_01_rrf_core.py`)

| Test | Input | Acceptance Criteria |
|------|-------|---------------------|
| `test_rates_cli_regression` | `example.nwk` + outgroup | CLI produces `_rates.csv` / `_rates.nwk`; key node rates match golden |
| `test_times_cli_r3f_compat` | `example.nwk` + outgroup | Produces `_times.csv` / `_timetree.nwk`; root age matches golden |
| `test_rates_times_cli_outputs` | Same as above | `rates-times` produces `_report.json` with both rates and times blocks |
| `test_times_normalize` | `--normalize` | Node ages are globally scaled while topology is preserved |
| `test_arithmetic_mean_cli` | `--mean arithmetic` | Runs successfully with arithmetic mean |
| `test_guard_options` | `--no-guard` / default | Root age differs when guard is disabled |
| `test_nexus_roundtrip` | NEXUS output re-read | NEXUS read/write preserves node count |

### 6.2 Calibration & CI (`test_02_calibration_ci.py`)

| Test | Input | Acceptance Criteria |
|------|-------|---------------------|
| `test_calibrate_bounds_cli` | Calibrations + `--method bounds` | Produces `_calibrated_*.csv`; time factor > 0 |
| `test_calibrate_effective_cli` | Calibrations + `--method effective` | Effective method runs successfully |
| `test_ci_default` | Calibrated result | CI file produced with lower/upper node-age bounds |
| `test_ci_with_n_sites` | `--n-sites` | CI width changes with n-sites |
| `test_ci_with_branch_var` | `--branch-var` | Runs using branch-length variances |
| `test_calibration_api_matches_cli` | Same input via API and CLI | Time factors are identical |
| `test_conflicting_calibration_rejected` | Conflicting bounds | CLI exits non-zero and reports conflict |
| `test_density_only_requires_effective` | Density-only calibration + effective | Runs successfully |
| `test_pinned_node_ci_flagged` | Pinned calibration | Corresponding node is flagged as pinned in CI table |

### 6.3 CorrTest / ddBD (`test_03_corrtest_ddbd.py`)

| Test | Input | Acceptance Criteria |
|------|-------|---------------------|
| `test_corrtest_cli_matches_golden` | `example.nwk` | CorrScore matches golden |
| `test_corrtest_api_full_precision` | Same | Full-precision API output matches golden |
| `test_corrtest_resample` | `--sister-resample 50` | Result matches `sr50` golden |
| `test_ddbd_cli_matches_golden` | Same | Birth/death rates match golden |
| `test_ddbd_anchor_matches_golden` | `--anchor-node` / `--anchor-time` | Anchor result matches golden |
| `test_ddbd_plot` | Same | Produces `_ddbd.png` |
| `test_ddbd_api_matches_cli` | Same input via API and CLI | Numerical outputs match |

### 6.4 tree2table / monophyly (`test_04_tree2table_monophyly.py`)

| Test | Input | Acceptance Criteria |
|------|-------|---------------------|
| `test_tree2table_branch_lengths` | `example.nwk` | CSV contains original branch lengths |
| `test_tree2table_node_ages` | RRF timetree | CSV node ages are non-empty |
| `test_monophyly_auto_groups` | `--taxa-table` | Auto-grouping results are as expected |
| `test_monophyly_taxon_table_groups` | Same | Monophyly verdicts per table group are correct |
| `test_monophyly_ad_hoc_tips` | `--tips` | Ad-hoc tip-set verdicts are correct |

### 6.5 CLI/API parity (`test_05_cli_api_parity.py`)

| Test | Input | Acceptance Criteria |
|------|-------|---------------------|
| `test_api_write_matches_cli_rates_times` | API-generated files vs CLI | DataFrame difference is zero |
| `test_report_json_merges_blocks` | Multi-stage analysis | `_report.json` merges all stage parameters |
| `test_cli_parameters_recorded` | Parameterized run | Report parameters match command-line flags |

### 6.6 Input validation (`test_06_input_validation.py`)

| Test | Input | Acceptance Criteria |
|------|-------|---------------------|
| `test_negative_branch_length_rejected` | Negative branch-length tree | Exits non-zero; stderr contains "negative" |
| `test_duplicate_tip_rejected` | Duplicate-tip tree | Exits non-zero; stderr contains "duplicate" |
| `test_polytomy_rejected_by_default` | Polytomous tree | Fails by default; passes with `--resolve random` |
| `test_too_few_ingroup_tips_rejected` | Outgroup contains ingroup tips | Exits non-zero |
| `test_incomplete_outgroup_rejected` | Incomplete outgroup | Fails by default; passes with `--outgroup-check warn` |
| `test_empty_calibration_row_rejected` | Empty calibration row | Exits non-zero; stderr contains "bound" |
| `test_density_only_with_bounds_method_rejected` | Density-only + bounds | Exits non-zero; stderr contains "effective" |
| `test_missing_input_file` | Non-existent input | Exits non-zero; stderr contains "not found" |

### 6.7 External tools (`test_07_blb_external.py`, `test_08_megacc_external.py`)

| Test | Input | Acceptance Criteria |
|------|-------|---------------------|
| `test_blb_full_pipeline` | `small.fasta` + IQ-TREE | Produces `blb_summary.csv` and `blb_replicates.csv` |
| `test_blb_iqtree_args` | `--iqtree-arg` | Custom arguments run without error |
| `test_megacc_with_calibrations` | `example.nwk` + calibrations + megacc | Produces `_megacc_times.csv` |
| `test_megacc_without_calibrations` | `example.nwk` + megacc | Produces `_megacc_times.csv` |

## 7. Skip Conditions

- **IQ-TREE tests**: skipped when `IQTREE` is unset and neither `iqtree2` nor `iqtree` is on PATH.
- **MEGA-CC tests**: skipped when `MEGACC` is unset and `megacc` is not on PATH.

## 8. Risks and Assumptions

1. External-tool version differences may cause tiny numeric drift; golden tests use a fixed tolerance.
2. IQ-TREE's `-nt AUTO` can be very slow on some machines, so BLB smoke tests fix `-nt 1 -m GTR+G` and use a minimal alignment.
3. MEGA-CC is not available on the main conda channels and is treated as an optional external dependency.
4. The desktop GUI is not covered here: it ships as the independent
   `OpenRelTime-Studio` product and is tested in that project. This suite
   exercises only the engine's CLI and Python API.

## 9. Regression and Maintenance

- Add a corresponding `test_*.py` under `e2e_tests/cases/` for each major new feature.
- Update e2e tests and golden data whenever CLI arguments or output formats change.
- Before each release, run:
  ```bash
  python -m pytest e2e_tests/ -v
  ```
