# OpenRelTime End-to-End Test Suite

**Language: English** · [中文](README-zh.md)

This directory contains the **full-flow, real-data end-to-end tests** for OpenRelTime. Instead of mocking internal functions, these tests run the actual `openreltime` analysis modules through both the command-line interface (CLI) and the Python API, and compare outputs against the project's built-in golden results.

The test files are shipped with the package so users can run them independently after installation to verify that the example analysis workflow can be fully reproduced in their local environment.

## Directory Layout

```
e2e_tests/
├── README.md                 # This file
├── README-zh.md              # 中文说明
├── conftest.py               # Shared fixtures, paths, and CLI runner
├── data/
│   ├── example.nwk           # Example phylogeny (Newick)
│   ├── example_taxa.tsv      # Example taxon grouping table
│   ├── example_calibrations.tsv  # Example calibration constraints
│   ├── example.fasta         # Larger synthetic alignment (kept as optional)
│   ├── small.fasta           # Tiny synthetic alignment for BLB smoke tests
│   ├── golden/               # Golden results (copied from data/golden/r3f/)
│   └── bad/                  # Directory for bad-input fixtures (generated on demand)
├── cases/
│   ├── test_01_rrf_core.py
│   ├── test_02_calibration_ci.py
│   ├── test_03_corrtest_ddbd.py
│   ├── test_04_tree2table_monophyly.py
│   ├── test_05_cli_api_parity.py
│   ├── test_06_input_validation.py
│   ├── test_07_blb_external.py
│   └── test_08_megacc_external.py
└── outputs/                  # Test outputs (ignored by .gitignore)
```

## Environment Requirements

1. Activate the `openreltime` conda environment:
   ```bash
   conda activate openreltime
   ```
2. Install the current package in editable mode with all optional dependencies:
   ```bash
   pip install --editable ".[dev,plot]"
   ```
3. (Optional) To run the external-tool tests, install the corresponding binaries:
   - **IQ-TREE** for the BLB tests. Install into the current environment with micromamba:
     ```bash
     micromamba install -p /path/to/openreltime/env -c conda-forge -c bioconda iqtree
     ```
   - **MEGA-CC `megacc`** for the megacc bridge tests. It is not available on the main conda channels; install manually and set the `MEGACC` environment variable.

## Test Data

| File | Source | Purpose |
|------|--------|---------|
| `data/example.nwk` | `data/examples/example.nwk` | Main phylogeny for RRF, calibrate, CI, CorrTest, ddBD, monophyly, etc. |
| `data/example_taxa.tsv` | `data/examples/example_taxa.tsv` | Taxon grouping for monophyly / tree2table tests |
| `data/example_calibrations.tsv` | `data/examples/example_calibrations.tsv` | Calibrations for calibrate, CI, and megacc tests |
| `data/small.fasta` | Randomly generated (10 taxa × 300 bp) | BLB / IQ-TREE external workflow smoke test |
| `data/example.fasta` | Randomly generated (274 taxa) | Optional larger alignment sample |
| `data/golden/*` | `data/golden/r3f/` | Regression baseline for CLI/API outputs |

All test data reside under `e2e_tests/data/`. Tests never create or modify files outside the repository root.

## Running the Tests

```bash
# Full suite (recommended)
python -m pytest e2e_tests/ -v

# RRF core tests only
python -m pytest e2e_tests/cases/test_01_rrf_core.py -v

# Generate an HTML report
python -m pytest e2e_tests/ -v --html=e2e_tests/outputs/report.html
```

> Note: the default `testpaths` in `pyproject.toml` only includes `tests/`, so
> you must explicitly pass `e2e_tests/` when running this suite.

## Test Modules and Coverage

| Module | Cases | Coverage |
|--------|-------|----------|
| `test_01_rrf_core.py` | 7 | `rates` / `times` CLI, golden regression, output files, normalization, arithmetic mean, guard options, NEXUS round-trip |
| `test_02_calibration_ci.py` | 9 | `calibrate` bounds/effective, CI default/branch-var/n-sites, API vs CLI parity, conflicting calibrations, density-only constraints, pinned-node flagging |
| `test_03_corrtest_ddbd.py` | 7 | CorrTest golden regression, full-precision API, resampling; ddBD golden, anchor, plot, API/CLI parity |
| `test_04_tree2table_monophyly.py` | 5 | `tree2table` branch lengths and node ages; monophyly auto/table/ad-hoc groups |
| `test_05_cli_api_parity.py` | 3 | API-written outputs match CLI, `_report.json` merging, parameter recording |
| `test_06_input_validation.py` | 8 | Negative branch lengths, duplicate tips, polytomies, incomplete outgroups, empty calibration rows, density-only vs method mismatch, missing input files |
| `test_07_blb_external.py` | 2 | Full BLB pipeline (IQ-TREE) and custom IQ-TREE argument forwarding |
| `test_08_megacc_external.py` | 2 | MEGA-CC bridge with and without calibrations |
| **Total** | **43** | 8 modules, engine CLI + API only |

## Skip Conditions

- `test_07_blb_external.py` is skipped automatically when the `IQTREE` environment variable is not set and `iqtree2` / `iqtree` is not on PATH.
- `test_08_megacc_external.py` is skipped automatically when the `MEGACC` environment variable is not set and `megacc` is not on PATH.

## Outputs

All generated files are written under `e2e_tests/outputs/`, which is ignored by git. Each test case writes to a separate subdirectory named `e2e_tests/outputs/<pytest_tmp_name>/`.

## Notes

- Tests invoke `python -m openreltime.cli` via `subprocess.run` to exercise the same CLI entry point users see.
- Tests that involve randomness use a fixed `--seed` so results are reproducible.
- This suite covers the engine only. The desktop GUI ships as the independent
  `OpenRelTime-Studio` product and has its own test suite there; no test here
  imports or launches it.
