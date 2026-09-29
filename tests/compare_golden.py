"""Compare OpenRelTime output against stored R3F golden CSVs (clade-matched).

Helper used by the pytest regression suite (tests/test_regression_golden.py) and by
``reproduce/`` scripts.  Run standalone:

    python tests/compare_golden.py <tree.nwk> <golden_dir> <outgroup,comma,sep>
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from openreltime.times import rrf_times  # noqa: E402
from openreltime.treeio import read_tree  # noqa: E402


def clade_map_from_csv(df: pd.DataFrame) -> dict[int, frozenset]:
    clades: dict[int, frozenset] = {}
    pending: list[int] = []
    for _, r in df.iterrows():
        nid = int(r["NodeId"])
        if r["NodeLabel"] != "-":
            clades[nid] = frozenset([r["NodeLabel"]])
        else:
            pending.append(nid)
    rows = {int(r["NodeId"]): r for _, r in df.iterrows()}
    while pending:
        progressed = False
        for nid in list(pending):
            row = rows[nid]
            try:
                d1, d2 = int(row["Des1"]), int(row["Des2"])
            except (TypeError, ValueError):
                continue
            if d1 in clades and d2 in clades:
                clades[nid] = clades[d1] | clades[d2]
                pending.remove(nid)
                progressed = True
        if not progressed:
            raise RuntimeError("cannot resolve golden clades from CSV")
    return clades


def golden_values(df: pd.DataFrame, column: str) -> dict[frozenset, float]:
    clades = clade_map_from_csv(df)
    out: dict[frozenset, float] = {}
    for _, r in df.iterrows():
        nid = int(r["NodeId"])
        if r["NodeLabel"] != "-":
            out[clades[nid]] = 0.0 if column == "Time" else float(r[column])
        else:
            out[clades[nid]] = float(r[column])
    return out


def compare(
    tree_path: Path,
    golden_dir: Path,
    outgroup: list[str],
    *,
    with_rate: bool = False,
) -> dict[str, float]:
    tree = read_tree(tree_path, outgroup=outgroup)
    res = rrf_times(tree)
    myclades = {n.node_id: n.clade for n in tree.walk()}

    csv = golden_dir / ("example_RRF_table.csv" if with_rate else "example_RRF_times.csv")
    gdf = pd.read_csv(csv)
    gtime = golden_values(gdf, "Time")
    mtime = {myclades[nid]: t for nid, t in res.times.items()}
    common = set(gtime) & set(mtime)
    gs = np.array([gtime[c] for c in common])
    ms = np.array([mtime[c] for c in common])
    mask = gs > 0
    slope = float(np.polyfit(gs[mask], ms[mask], 1)[0])
    corr = float(np.corrcoef(gs[mask], ms[mask])[0, 1])
    rel = float(np.max(np.abs(gs[mask] - ms[mask]) / gs[mask]))
    median_rel = float(np.median(np.abs(gs[mask] - ms[mask]) / gs[mask]))
    out = {
        "time_slope": slope,
        "time_corr": corr,
        "time_max_rel": rel,
        "time_median_rel": median_rel,
        "n_common": len(common),
    }

    if with_rate:
        rate_csv = golden_dir / "example_RRF_rates.csv"
        rdf = pd.read_csv(rate_csv)
        grate = golden_values(rdf, "Rate")
        mrate = {myclades[nid]: t for nid, t in res.rates.items()}
        rcommon = set(grate) & set(mrate)
        gs2 = np.array([grate[c] for c in rcommon])
        ms2 = np.array([mrate[c] for c in rcommon])
        mask2 = gs2 > 0
        out["rate_slope"] = float(np.polyfit(gs2[mask2], ms2[mask2], 1)[0])
        out["rate_median_rel"] = float(
            np.median(np.abs(gs2[mask2] - ms2[mask2]) / gs2[mask2])
        )
        out["rate_max_rel"] = float(
            np.max(np.abs(gs2[mask2] - ms2[mask2]) / gs2[mask2])
        )
    return out


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print(__doc__)
        return 2
    tree_path = Path(argv[1])
    golden_dir = Path(argv[2])
    outgroup = argv[3].split(",")
    result = compare(tree_path, golden_dir, outgroup, with_rate=True)
    for key, value in result.items():
        print(f"{key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
