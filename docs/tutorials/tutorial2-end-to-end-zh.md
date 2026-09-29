# 教程 2：带校正与置信区间的端到端定年

**语言：中文** · [English](tutorial2-end-to-end.md)

流程：比对 → IQ-TREE → OpenRelTime → 校正 → 置信区间。

## 校正与解析置信区间

```bash
# 1. 由比对推断系统树（IQ-TREE）
iqtree -s alignment.fasta -m GTR+G -B 1000 -pre mam
# 带枝长的树是 mam.treefile
```

```python
import openreltime as ort
from openreltime.calibrate import parse_calibrations
from openreltime.ci import confidence_interval

tree = ort.read_tree("mam.treefile", outgroup=["Out1", "Out2"])
times = ort.rrf_times(tree)

# 2. 从 TSV 文件读取校正（按末端集合的 MRCA 或内部节点 id）
cals = parse_calibrations("calibrations.tsv")

# 3a. 硬边界 → 绝对时间
cal = ort.calibrate(times, cals, method="bounds")
cal.write("mam_bounds")

# 3b. 密度 → 有效边界（msz236）+ 经验节点年龄区间
cal_eff = ort.calibrate(times, cals, method="effective",
                        n_effective=10_000, seed=42)
cal_eff.write("mam_effective")   # 写出 mam_effective_effective_bounds.csv，
                                 # 其 2.5/97.5 两列即经验置信区间

# 4. 解析置信区间（delta method，msz236）——此处针对 bounds 那次结果
ci = confidence_interval(cal, seq_length=1000)
ci.write("mam_ci")
```

`calibrations.tsv` 示例（制表符分隔；`.` 表示空单元格，每行都要含齐六个字段）：

```
node_id	taxon_set	min_bound	max_bound	density	density_params
-	Homo_sapiens|Pan_troglodytes	6.0	8.5	.	.
-	Elephas_maximus|Loxodonta_africana	.	.	exponential	offset=60;mean=20
-	Dasypus_novemcinctus|Choloepus_didactylus	70.0	105.0	.	.
```

这些行只演示文件格式，边界数值是随意填写的示意值，并非建议的化石约束；见
`data/examples/README-zh.md`。

> 第 4 步求导的对象，是传给 `confidence_interval` 的那个结果，这里传入的是
> **bounds** 的结果。`confidence_interval` 会复现所接收结果的估计器——平均约定、
> 速率比保护以及*实际使用*的边界——因此把 `cal_eff` 交给它同样成立：effective
> 结果会按其有效边界求导，不会退回原始的（可能只有密度的）那几行。
>
> 两条路径回答的问题不同。解析 CI 把枝长噪声传播为节点年龄的不确定度；有效边界
> 方法自己的 2.5 与 97.5 经验分位数（`mam_effective_effective_bounds.csv`）刻画
> 的是校正密度本身。两者不一致时，建议同时报告。另外注意：`mam_ci.csv` 会给钉在
> 硬边界上的节点标记 `se_reliable = False`，因为该处式 (7) 的导数恰为 0。

## BLB：计入树不确定性的置信区间（需要 IQ-TREE）

```bash
openreltime blb -a alignment.fasta --iqtree iqtree \
    --outgroup "Out1,Out2" --gamma 0.7 \
    --n-subsample 10 --n-replicate 10 \
    --iqtree-arg=-m --iqtree-arg=GTR+G --seed 42 -o blb_out
```

`blb_out/blb_summary.csv` 逐 clade 给出中位时间、2.5 与 97.5 分位数，以及包含该
clade 的复制比例（`topology_ok`）。`blb_replicates.csv` 保留每一次复制的完整结果
表。

上一步：[教程 1](tutorial1-quickstart-zh.md)；下一步：
[教程 3](tutorial3-priors-megacc-zh.md)。
