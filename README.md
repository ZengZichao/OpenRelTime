# OpenRelTime

**Language: English** · [中文](README-zh.md)

**OpenRelTime** is an open-source Python implementation of the **relative
rate framework (RRF)** for non-Bayesian molecular dating (Tamura et al.
2018, *MBE* 35:1770-1782), the method behind RelTime in MEGA.  It estimates
relative lineage rates and divergence times from a branch-length tree,
converts calibrations (min/max bounds or probability densities) into
absolute times with analytical confidence intervals, and adds CorrTest
rate-autocorrelation testing, birth-death (ddBD) tree-prior estimation and a
bag-of-little-bootstraps (BLB) pipeline — everything scriptable from Python
or the command line.

- Documentation (English): [docs/README.md](docs/README.md) — user guide,
  parameter handbook, methods, validation, tutorials
- 中文文档（全手册亦有中文版）：[README-zh.md](README-zh.md) ·
  [文档目录](docs/README.md) · [使用说明](docs/usage-zh.md)
- API/CLI equivalence, reproducible golden standards, GPL-3.0

```
@software{openreltime,
  author = {曾, 子超},
  title  = {OpenRelTime: relative rate framework based molecular dating in Python},
  year   = {2026},
  url    = {https://github.com/ZengZichao/OpenRelTime},
  note   = {v0.1.0}
}
```

## Why OpenRelTime

| Tool | Language | Open source | Calibration | Analytical CI | CorrTest | ddBD | API + CLI |
|---|---|---|---|---|---|---|---|
| MEGA11/12 (RelTime) | Pascal | GPL-3 | ✓ | ✓ | ✓ | – | weak |
| R3F | R | GPL-3 | – | – | ✓ | ✓ | – |
| **OpenRelTime** | **Python** | **GPL-3** | **✓ (bounds + densities)** | **✓** | **✓** | **✓** | **✓** |

OpenRelTime reproduces R3F numerically on the 274-tip mammalian benchmark
(regression slope 1.0000; observed median relative deviation 6.8e-16 for rates
and 8.7e-16 for times).  The automated gate G2a/G2b asserts a slope >= 0.999
and a median relative deviation <= 1e-6; the largest observed deviation
(~4.8e-15 rates, 3.6e-14 times) is measured but not asserted by the gate — see
`docs/validation.md`.  On top of the shared core OpenRelTime adds the modules
R3F lacks — calibration conversion and analytical confidence intervals (msz236)
— plus an end-to-end BLB pipeline.

## Installation

> **Distribution.** This release ships as source plus the sdist and wheel
> attached to [the v0.1.0 release](https://github.com/ZengZichao/OpenRelTime/releases/tag/v0.1.0).
> It is not on a package index yet; an index release will be announced here.

From a checkout:

```bash
git clone https://github.com/ZengZichao/OpenRelTime
cd OpenRelTime
pip install --editable ".[dev,plot]"   # core (numpy, scipy, pandas, click)
                                       # + matplotlib visualisation ([plot])
                                       # + pytest/ruff/mypy ([dev], contributors)
```

Or install the released wheel directly from GitHub — no clone, and the `[plot]`
extra brings matplotlib:

```bash
pip install "openreltime[plot] @ https://github.com/ZengZichao/OpenRelTime/releases/download/v0.1.0/openreltime-0.1.0-py3-none-any.whl"
```

The desktop GUI is a separate product, installed separately (see below).

`plot` and `dev` are the only extras the package declares; combine them in one
call, e.g. `".[plot,dev]"`.

Requirements: Python >= 3.10, Linux/macOS/Windows.  There are no compiled
extensions; installation takes seconds.  CI runs the test suite on Linux and
macOS across Python 3.10-3.13 (`.github/workflows/ci.yml`); Windows is a
declared target but is not exercised by CI (see `docs/validation.md`).

## OpenRelTime Studio (desktop GUI) — separate product

**OpenRelTime Studio** is a native desktop app for researchers who prefer
point-and-click: load a Newick/NEXUS tree, run RRF, add calibrations by clicking
internal nodes on the canvas, then run calibration, analytical CIs, CorrTest and
ddBD — all through this engine's public Python API, with CSV/NEXUS/JSON/PNG
export and a reproducible equivalent-command-line panel.
Its interface is fully bilingual (English / 简体中文) and switches live, it ships
light and dark themes, a hand-drawn vector SVG icon, and a built-in example tree
with demo calibrations so first-time users need no data.

It is **not** part of this repository and is **not** installed from it.  It ships
as its own project and depends on the engine as an ordinary third-party package.
Neither is on an index yet, so install both wheels from their releases:

```bash
pip install \
  "openreltime[plot] @ https://github.com/ZengZichao/OpenRelTime/releases/download/v0.1.0/openreltime-0.1.0-py3-none-any.whl" \
  "OpenRelTime-Studio @ https://github.com/ZengZichao/OpenRelTime-Studio/releases/download/v0.1.0/openreltime_studio-0.1.0-py3-none-any.whl"
openreltime-studio
```

On Apple Silicon, Studio's release also attaches a prebuilt
`OpenRelTimeStudio-v0.1.0-macOS-arm64.zip` that needs no Python at all.

Home page: <https://github.com/ZengZichao/OpenRelTime-Studio> ·
GUI manual: see that project's documentation.  Installing the engine alone leaves
the Python API and every CLI subcommand above unchanged — Studio adds a front
end, never new analysis code.

## Quick start (30 seconds)

```python
import openreltime as ort

tree = ort.read_tree("example.nwk", outgroup=["Out1", "Out2"])
times = ort.rrf_rates_times(tree)             # rates + relative times
times.write("quickstart", with_rate=True, nexus=True)   # CSV + NEXUS timetree
```

Equivalent CLI:

```bash
openreltime rates-times -i example.nwk --outgroup "Out1,Out2" -o quickstart
```

From an alignment (BLB, requires IQ-TREE):

```bash
openreltime blb -a alignment.fasta --iqtree iqtree --seed 42 -o blb_out
```

## Features

- **RRF core** — geometric mean (R3F/msy044-compatible) and arithmetic mean
  (PNAS 2012) conventions; extreme-rate guard; polytomy handling; Newick +
  NEXUS (`[&rate=...]`) I/O.
- **Calibration** — hard min/max bounds and probability densities
  (uniform/exponential/normal/lognormal) converted to absolute times via the
  global factor *f*; the effective-bounds resampling method of msz236.
- **Analytical CIs** — delta method (msz236 eqs 7-14) with flexible sampling
  variance sources and truncation at hard bounds; empirical node-age CIs
  from the effective-bounds replicates.
- **CorrTest** — autocorrelated-rates test with the published
  fixed-coefficient logistic model.
- **ddBD** — birth-death speciation prior (birth, death, sampling fraction)
  for Bayesian dating pipelines.
- **BLB pipeline** — two-level site resampling through IQ-TREE with
  per-clade time aggregation.
- **MEGA-CC bridge** (optional) — generate `.mao` files, run `megacc`, parse
  RelTime output back to OpenRelTime tables.

## Command-line overview

```
openreltime rates        -i tree.nwk [--outgroup og.txt] [--mean geometric] [-o prefix]
openreltime times        -i tree.nwk [--normalize] [--no-guard] [--r3f-compat] [-o prefix]
openreltime rates-times  -i tree.nwk [--plot timetree.png] [--r3f-compat] [-o prefix]
openreltime calibrate    -i tree.nwk -c calibrations.tsv [--method effective] [-o prefix]
openreltime ci           -c <calibrated prefix> [--n-sites 1000] [-o prefix]
openreltime corrtest     -i tree.nwk [--sister-resample 100] [-o prefix]
openreltime ddbd         -i tree.nwk [--anchor-time 1.85] [-o prefix]
openreltime tree2table   -i tree.nwk [--time] [-o prefix]
openreltime monophyly   -i tree.nwk [--taxon-table taxa.tsv] [--groups A,B] [-o prefix]
openreltime blb          -a aln.fasta --iqtree iqtree [--gamma 0.7 ...] [-o dir]
openreltime megacc       -i tree.nwk -c calibrations.tsv --megacc /path/to/megacc
```

Every analysis stage whose parameters can move the numbers (`rates`, `times`,
`rates-times`, `calibrate`, `ci`, `corrtest`, `ddbd`) merges its parameters,
warnings and random seed into `<prefix>_report.json`; runs that share a prefix
contribute their own key. `tree2table`, `monophyly`, `blb` and `megacc` write
only their table.

## Calibration file (`calibrations.tsv`)

```
node_id	taxon_set	min_bound	max_bound	density	density_params
	Homo_sapiens|Pan_troglodytes	6.0	8.5	.
	Elephas_maximus|Loxodonta_africana	.	.	exponential	offset=60;mean=20
	Dasypus_novemcinctus|Choloepus_didactylus	70	105	.
```

Identify nodes by MRCA of a taxon set *or* by internal node id; give hard
bounds, a density, or both.  Densities: `uniform`, `exponential`,
`normal`, `lognormal` (parameters as `key=value;...`).

## Validation

The test suite encodes the project's acceptance gates; each gate, what it
actually asserts, and which shipped modules it does *not* cover are listed in
[`docs/validation.md`](docs/validation.md):

```bash
python -m pytest tests/
```

- **G1** hand-computed 3-/4-taxon RRF values (both conventions, 1e-12)
- **G2a/G2b** R3F regression on the 274-tip mammalian tree (gate asserts slope
  >= 0.999 and median relative deviation <= 1e-6; observed medians 6.8e-16
  rates / 8.7e-16 times)
- **G2c** NEXUS `[&rate=...]` interop round-trip
- **G6/G7** CorrTest agreement with the stored R3F CorrScore (0.9996);
  calibration bounds/conflict handling
- Golden files + example inputs + the R script that produced the goldens ship
  in `data/golden/r3f/`, `data/examples/` (see its `README.md`),
  `data/PROVENANCE.md` and `reproduce/golden_r3f.R`.

## Documentation

The full manual set exists in **both English and Chinese**; each document links
to its counterpart at the top.  Index: [docs/README.md](docs/README.md).

| Document | Content | 中文 |
|---|---|---|
| [docs/usage-en.md](docs/usage-en.md) | Detailed user guide — every subcommand, output file, FAQ and error message | [docs/usage-zh.md](docs/usage-zh.md) |
| [docs/parameters.md](docs/parameters.md) | Parameter handbook — every option, with a per-subcommand option matrix | [docs/parameters-zh.md](docs/parameters-zh.md) |
| [docs/methods.md](docs/methods.md) | Methods, assumptions and limitations | [docs/methods-zh.md](docs/methods-zh.md) |
| [docs/validation.md](docs/validation.md) | Validation gates, coverage gaps, tested platforms | [docs/validation-zh.md](docs/validation-zh.md) |
| [docs/tutorials/](docs/tutorials/) | Three step-by-step tutorials | [same directory](docs/tutorials/), `*_zh.md` |
| [docs/adr/README.md](docs/adr/README.md) | Architecture-decision-record index | [docs/adr/README-zh.md](docs/adr/README-zh.md) |
| [CHANGELOG.md](CHANGELOG.md) | Release notes | [CHANGELOG-zh.md](CHANGELOG-zh.md) |
| [THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md) | Reference material and adapted constants, per source | [THIRD-PARTY-NOTICES-zh.md](THIRD-PARTY-NOTICES-zh.md) |
| [data/PROVENANCE.md](data/PROVENANCE.md) | Origin, licence and citation of every bundled third-party file | [data/PROVENANCE-zh.md](data/PROVENANCE-zh.md) |

## License and reuse

GPL-3.0-or-later (`LICENSE`).  Behavioural references (R3F, MEGA sources,
ape) are credited in `THIRD-PARTY-NOTICES.md` with per-file provenance in
`data/PROVENANCE.md`.  The CorrTest repository (no licence) was not copied;
only published constants of Tao et al. (2019) are used.

## Authors

曾子超 (Zichao Zeng) · [zengzichao@sjtu.edu.cn](mailto:zengzichao@sjtu.edu.cn) ·
[ORCID 0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)

The preferred citation is recorded in [`CITATION.cff`](CITATION.cff).
