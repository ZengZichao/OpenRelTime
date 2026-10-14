# Example data

**Language: English** · [中文](README-zh.md)

Files in this directory are **worked examples** used by the tutorials and the
test suite.  They are input fixtures, not analysis results, and they are not
recommended settings for a real study.

## `example.nwk`

A 274-tip mammalian chronogram (7,370 sites) redistributed from R3F
(Tao, Sharma, Tamura & Kumar 2025), originally dos Reis et al. 2012.  Used as
the regression input for the golden-standard tests.  Provenance and licence:
see `data/PROVENANCE.md`.  The outgroup used by the tests and
tutorials is the three monotreme tips `Ornithorhynchus_anatinus`,
`Zaglossus_bruijni`, `Tachyglossus_aculeatus`.

## `example_calibrations.tsv` and `example_taxa.tsv`

**These files illustrate the file format only.**  They show how a
`key=value` calibration table and a tip->group classification table are
laid out; the *numbers* in them must not be read as suggested fossil
constraints or as a curated taxonomy.

**`example_calibrations_feasible.tsv` is the runnable counterpart**: four
mutually consistent hard bounds that solve under the documented default
`method='bounds'` (the file above cannot, because its second row is
density-only).

In particular the three bounds in `example_calibrations.tsv`

| taxon set | bound used here |
|---|---|
| `Homo_sapiens|Pan_troglodytes` | min 6.0 / max 8.5 |
| `Elephas_maximus|Loxodonta_africana` | exponential density `offset=60;mean=20` |
| `Dasypus_novemcinctus|Choloepus_didactylus` | min 70.0 / max 105.0 |

are arbitrary placeholders chosen to exercise each parser branch (hard bounds,
a density, and a wide hard bound).  Taken literally they are mutually
inconsistent: the `Elephas`/`Loxodonta` density implies a much older
`Dasypus`/`Choloepus` node than the 70-105 bound allows, so the solver must
rescale lineage rates by roughly a factor of 23 to make the file feasible.
That infeasibility is a property of the *placeholder values*, not of the
format; supply real, individually justified bounds for any analysis you intend
to report.

**What running this file actually does.**  Under the documented default
`method='bounds'` it is refused before the solver starts, because the second
row is density-only and the parser names that row.  Under
`method='effective'` it is refused *after* the solver converges.  Measured on
the 274-tip benchmark at the msz236 protocol (`n_effective=10000`,
`seed=20260906`), the effective-bounds solve lands on

| quantity | value |
|---|---|
| `time_factor` | 4432.82 |
| root age (node 272) | 1162.86 Ma |
| youngest node (node 307) | **-3998.56 Ma** |
| parent-child inversions | 4 (node 306 at 1746.23 Ma below its ancestor 305 at 105 Ma; tips 32 and 33 below node 307; node 372 at 75.49 Ma below its ancestor 371 at 8.5 Ma) |
| `ImpliedRate` entries that are `NaN` | 4 |

A negative age and a descendant older than its ancestor are not alternative
answers; they are the absence of one.  `openreltime.calibrate` now checks the
converged vector against the tree and **raises** `ValueError`, naming the
offending node ids, instead of returning those ages with four log warnings as
0.1.0 did (`calibrate.py`'s repair loop tested only the worst violator of each
round, so an inversion created elsewhere in the tree was never re-examined).
The exact figures above are what the check rejects, and they vary with the seed
— which is itself evidence that the cause is structural rather than a one-seed
artefact.  Treat this file as a **format illustration that is expected to fail**,
and use `example_calibrations_feasible.tsv` below for anything you mean to run.

> Why the warning lives in this README rather than in the TSV header: the
> loader (`openreltime.calibrate.parse_calibrations`) reads the file with a
> strict tab-separated `key=value` parser and has **no comment syntax**.  A
> leading `#` line is parsed as a data row and makes the whole file fail to
> load (it is read as the header's first field, after which the required
> `min_bound`/`max_bound` columns are missing).  A comment header is therefore
> not a safe place for this notice, and `data/examples/README.md` is used
> instead.

## `example_calibrations_feasible.tsv` — the one file that runs as shipped

Unlike the format-illustration table above, this file is a **worked example**:
four hard bounds, no densities, chosen so that all four are jointly satisfiable
by a single global time factor.  It therefore runs under the documented default
`method='bounds'`, which `example_calibrations.tsv` cannot (its density-only
second row is rejected before the solver starts):

```bash
openreltime calibrate -i data/examples/example.nwk \
    -c data/examples/example_calibrations_feasible.tsv \
    --outgroup "Ornithorhynchus_anatinus,Zaglossus_bruijni,Tachyglossus_aculeatus"
```

```python
cal = openreltime.calibrate(times, openreltime.parse_calibrations(
    "data/examples/example_calibrations_feasible.tsv"))          # method='bounds'
```

Both give `time_factor = 371.60` and these node ages, each strictly inside its
own bound, with no warning in the result and no missing implied rate:

| taxon set | bound | calibrated age |
|---|---|---|
| `Homo_sapiens\|Pan_troglodytes` | 5.0 – 8.0 | 5.59 |
| `Loxodonta_africana\|Elephas_maximus` | 2.0 – 4.0 | 2.70 |
| `Mus_musculus\|Rattus_norvegicus` | 12.0 – 18.0 | 14.76 |
| `Dasypus_novemcinctus\|Choloepus_didactylus` | 43.0 – 52.0 | 47.50 |

The bounds are still illustrative soft bounds, not a curated fossil
calibration set; they were chosen for consistency with each other, not for
palaeontological justification.  Three of the four taxon sets are two-tip sets
whose MRCA also subtends close relatives (`Pan_paniscus`; the other
`Mus`/`Rattus` species; the other xenarthrans), so the loader prints the
warn-only "not a complete clade" notice while still targeting the correct MRCA
— the standard behaviour of an MRCA calibration, not an error.
