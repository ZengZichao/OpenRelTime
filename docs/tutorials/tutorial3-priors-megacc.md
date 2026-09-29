# Tutorial 3: MCMCTree priors with ddBD and the MEGA-CC bridge

**Language: English** · [中文](tutorial3-priors-megacc-zh.md)

## Birth-death prior for Bayesian dating

```python
import openreltime as ort

tree = ort.read_tree("mam.treefile", outgroup=["Out1", "Out2"])

# sampling fraction fixed (only birth/death estimated)
bd = ort.ddbd(tree, sampling_frac=0.5, anchor_time=1.85)

# everything estimated (default)
bd_full = ort.ddbd(tree, anchor_time=1.85, measure="KL")

print(bd.birth_rate, bd.death_rate, bd.sampling_frac)
bd.write("mam_ddbd")
```

The estimated birth rate, death rate and sampling fraction parameterise a
birth-death node-age density usable as an MCMCTree/BEAST tree prior.  Keep
`anchor_time` (and therefore the maximum node age) <= 10 time units.

Note the semantics of `sampling_frac`, which follow R3F: it sets the fraction
*reported*, while the likelihood fit always leaves rho free, so the value does
not constrain the fit.  The free-rho optimum stays visible as
`sampling_frac_fitted` in the report.

## MEGA-CC bridge (optional extra)

For migration from MEGA workflows or third-party cross-checks:

```python
from openreltime.megacc import run_megacc, parse_megacc_output

mao, outputs = run_megacc(
    "mam.treefile",
    megacc="/opt/megacc",
    calibrations="calibrations.tsv",
    outgroup=["Out1", "Out2"],
    workdir="megacc_run",
)
table = parse_megacc_output(outputs, outgroup=["Out1", "Out2"])
```

`run_megacc` generates a minimal `reltimeFromBranchLengths.mao`
(`ppRelTimeBLens = true` plus optional overrides), probes the megacc version
banner and validates every argument before starting the process (no shell).
The outputs are parsed back into an R3F-compatible node-time table so MEGA
results can be compared node-by-node with `openreltime.rrf_times`.

The equivalent command line:

```bash
openreltime megacc -i mam.treefile -c calibrations.tsv --megacc /opt/megacc \
    --outgroup "Out1,Out2" --workdir megacc_run -o mam_mega
# results land in mam_mega_megacc_times.csv
```

Two limits apply (see `docs/validation.md`): the bridge needs a real `megacc`
executable, which CI neither installs nor runs, so MEGA-CC agreement is
established against the published R3F/MEGA results and the shipped golden files
rather than by a head-to-head run here; and OpenRelTime's analytical RRF layer
lacks MEGA's likelihood-ratio refinement, so node-by-node differences are
expected — this is a comparison tool, not an equivalence claim.

Previous: [Tutorial 2](tutorial2-end-to-end.md).
