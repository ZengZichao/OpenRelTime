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
carries all six fields):

```
node_id	taxon_set	min_bound	max_bound	density	density_params
-	Homo_sapiens|Pan_troglodytes	6.0	8.5	.	.
-	Elephas_maximus|Loxodonta_africana	.	.	exponential	offset=60;mean=20
-	Dasypus_novemcinctus|Choloepus_didactylus	70.0	105.0	.	.
```

These rows illustrate the file format; the bounds are arbitrary values, not
suggested fossil constraints.  See `data/examples/README.md`.

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
