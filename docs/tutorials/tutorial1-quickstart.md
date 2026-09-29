# Tutorial 1: five minutes with OpenRelTime

**Language: English** · [中文](tutorial1-quickstart-zh.md)

```python
import openreltime as ort

# 1. read a branch-length tree; root on an outgroup and prune it
tree = ort.read_tree(
    "data/examples/example.nwk",
    outgroup=[
        "Ornithorhynchus_anatinus",
        "Zaglossus_bruijni",
        "Tachyglossus_aculeatus",
    ],
)

# 2. relative lineage rates and node times
rates = ort.rrf_rates(tree)
times = ort.rrf_rates_times(tree)          # both, one pass

# 3. write CSV tables, a timetree and a rate-annotated NEXUS file
times.write("tutorial1", with_rate=True, nexus=True)

# 4. is the clock autocorrelated?
corr = ort.corrtest(tree)
print(corr.score, corr.p_band)

# 5. birth-death speciation prior from the relative times
bd = ort.ddbd(tree, anchor_time=1.85)
print(bd.birth_rate, bd.death_rate, bd.sampling_frac)
```

The same analyses from a shell:

```bash
openreltime rates-times -i example.nwk --outgroup "Ornithorhynchus_anatinus,Zaglossus_bruijni,Tachyglossus_aculeatus" -o tutorial1 --plot tutorial1.png
openreltime corrtest -i example.nwk --outgroup "Ornithorhynchus_anatinus,Zaglossus_bruijni,Tachyglossus_aculeatus" -o tutorial1
openreltime ddbd -i example.nwk --outgroup "Ornithorhynchus_anatinus,Zaglossus_bruijni,Tachyglossus_aculeatus" --anchor-time 1.85 -o tutorial1
```

Files written (all share the `tutorial1` prefix; the three commands merge their
parameters, warnings and seed into the single `tutorial1_report.json`):

| File | From |
|---|---|
| `tutorial1_rates_times.csv` | `rates-times` (rate + time columns; `--r3f-compat` would rename it `_RRF_table.csv`) |
| `tutorial1_timetree.nexus` | `rates-times` (relative-time tree with `[&rate=...]` annotations, opens in FigTree) |
| `tutorial1.png` | `rates-times --plot` (rate-coloured timetree image) |
| `tutorial1_corrtest.txt` | `corrtest` (CorrScore + P band) |
| `tutorial1_ddbd.txt` | `ddbd` (birth / death / sampling rates) |
| `tutorial1_report.json` | every run (parameters + warnings + seed, for reproducibility) |

Note that `rates-times` writes a `_rates_times.csv` table and a `_timetree.nexus`
tree; it does **not** emit a separate `_rates.csv` or a `_timetree.nwk`.  For a
standalone rates table use the `rates` subcommand (`<prefix>_rates.csv` +
`<prefix>_rates.nwk`), and for a plain Newick timetree use `times`
(`<prefix>_times.csv` + `<prefix>_timetree.nwk`).

Next: [Tutorial 2](tutorial2-end-to-end.md).
