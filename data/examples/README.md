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

> Why the warning lives in this README rather than in the TSV header: the
> loader (`openreltime.calibrate.parse_calibrations`) reads the file with a
> strict tab-separated `key=value` parser and has **no comment syntax**.  A
> leading `#` line is parsed as a data row and makes the whole file fail to
> load (it is read as the header's first field, after which the required
> `min_bound`/`max_bound` columns are missing).  A comment header is therefore
> not a safe place for this notice, and `data/examples/README.md` is used
> instead.
