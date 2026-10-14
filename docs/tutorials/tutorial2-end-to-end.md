# Tutorial 2: end-to-end dating with calibrations and CIs

**Language: English** · [中文](tutorial2-end-to-end-zh.md)

Workflow: alignment -> IQ-TREE -> OpenRelTime -> calibrations -> CI.

## Calibration and analytical CIs

```bash
# 1. infer a phylogeny from an alignment (IQ-TREE)
iqtree -s alignment.fasta -m GTR+G -B 1000 -pre mam
# the branch-length tree is mam.treefile
```

```python
import openreltime as ort
from openreltime.calibrate import parse_calibrations
from openreltime.ci import confidence_interval

tree = ort.read_tree("mam.treefile", outgroup=["Out1", "Out2"])
times = ort.rrf_times(tree)

# 2. calibrations from a TSV file (node MRCA or internal node id)
cals = parse_calibrations("calibrations.tsv")

# 3a. hard bounds -> absolute times
cal = ort.calibrate(times, cals, method="bounds")
cal.write("mam_bounds")

# 3b. densities -> effective bounds (msz236) + empirical node-age CIs
cal_eff = ort.calibrate(times, cals, method="effective",
                        n_effective=10_000, seed=42)
cal_eff.write("mam_effective")   # writes mam_effective_effective_bounds.csv,
                                 # whose 2.5/97.5 columns are the empirical CIs

# 4. analytical confidence intervals (delta method, msz236) for the bounds run
ci = confidence_interval(cal, seq_length=1000)
ci.write("mam_ci")
```

`calibrations.tsv` example (tab-separated; `.` marks an empty cell, every row
carries all six fields).  These four rows are a **worked set**, not just a
format sketch: three hard bounds plus one row carrying a density *and* its own
hard bounds, chosen so that a single global time factor can satisfy all four.
On the shipped 274-tip benchmark (`data/examples/example.nwk`, rooted on the
three monotremes) they give `f = 371.60` under `method="bounds"` and solve under
`method="effective"` too, with no warning on either run.  The bounds are still
illustrative soft bounds, not suggested fossil constraints.

```
node_id	taxon_set	min_bound	max_bound	density	density_params
-	Homo_sapiens|Pan_troglodytes	5.0	8.0	.	.
-	Mus_musculus|Rattus_norvegicus	12.0	18.0	.	.
-	Dasypus_novemcinctus|Choloepus_didactylus	43.0	52.0	.	.
-	Elephas_maximus|Loxodonta_africana	1.5	4.5	normal	mean=2.7;sd=0.3
```

Two rules decide whether a table like this runs at all, and both stop the solve
instead of quietly working around the problem:

* a row with a density but **no** hard bound constrains nothing under
  `method="bounds"` and is rejected when the table is parsed.  The fourth row
  therefore carries `1.5`/`4.5` as well as the density: step 3a uses those
  bounds, step 3b replaces them with the dates it samples;
* rows that no rate rescaling of the tree can satisfy **jointly** raise
  `ValueError`, naming the negative ages and the parent-child inversions.  In
  0.1.0 the same set came back as a result carrying warnings.  The previous,
  placeholder version of this table — a density-only `Elephas`/`Loxodonta` row
  beside a `70`–`105` `Dasypus`/`Choloepus` bound — was one such set, and the
  shipped `data/examples/example_calibrations.tsv` still is one: it is
  documented as a format illustration that is expected to fail, see
  `data/examples/README.md`.  For a table that runs as it stands use
  `data/examples/example_calibrations_feasible.tsv` — the same four hard bounds
  without the density row.

> Step 4 qualifies a **bounds** calibration, because that is the result it is
> handed.  `confidence_interval` reproduces the estimator behind whatever
> result it receives — its averaging convention, its guard and the bounds
> actually used — so feeding it `cal_eff` is also valid: an effective-bounds
> result is differentiated along its effective bounds rather than along the
> raw, possibly density-only rows.  The two routes answer different questions.
> The analytic CI propagates branch-length noise into node-age uncertainty;
> the effective method's own 2.5/97.5 empirical quantiles in
> `mam_effective_effective_bounds.csv` describe the calibration density itself.
> Report both when they differ, and note that a node pinned at a hard bound is
> flagged `se_reliable = False` in `mam_ci.csv` because its eq-(7) derivative
> is exactly zero there.

## BLB: tree-uncertainty confidence intervals (needs IQ-TREE)

```bash
openreltime blb -a alignment.fasta --iqtree iqtree \
    --outgroup "Out1,Out2" --gamma 0.7 \
    --n-subsample 10 --n-replicate 10 \
    --iqtree-arg=-m --iqtree-arg=GTR+G --seed 42 -o blb_out
```

`blb_out/blb_summary.csv` reports per-clade median times and 2.5/97.5
percentiles together with the fraction of replicates containing the clade
(`topology_ok`), and `blb_replicates.csv` keeps every replicate table.

Previous: [Tutorial 1](tutorial1-quickstart.md); next:
[Tutorial 3](tutorial3-priors-megacc.md).
