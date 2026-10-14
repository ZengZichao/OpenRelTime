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

`calibrations.tsv` 示例（制表符分隔；`.` 表示空单元格，每行都要含齐六个字段）。
下面四行是**一套真正能算的示例**，不只是格式示意：三条硬边界，外加一行同时给出
密度与自己的硬边界，四行选取的数值可以让同一个全局时间因子同时满足。在随包的
274 末端基准树上（`data/examples/example.nwk`，以三个单孔类末端为外群），它们在
`method="bounds"` 下给出 `f = 371.60`，在 `method="effective"` 下同样可解，两次
运行都没有警告。这些边界依旧是示意的软边界，不是建议采用的化石约束。

```
node_id	taxon_set	min_bound	max_bound	density	density_params
-	Homo_sapiens|Pan_troglodytes	5.0	8.0	.	.
-	Mus_musculus|Rattus_norvegicus	12.0	18.0	.	.
-	Dasypus_novemcinctus|Choloepus_didactylus	43.0	52.0	.	.
-	Elephas_maximus|Loxodonta_africana	1.5	4.5	normal	mean=2.7;sd=0.3
```

决定一张这样的表能不能跑的，是两条规则；两条都是在求解时停下来，而不是悄悄绕过去：

* 只给密度、**不给**硬边界的行，在 `method="bounds"` 下不构成任何约束，解析校正表
  时就会被拒绝。所以第四行除了密度还写了 `1.5`/`4.5`：第 3a 步用这两个边界，第 3b
  步把它们换成自己抽样出来的两个日期；
* 若一棵树无论怎样缩放速率都无法**同时**满足各行，就抛出 `ValueError`，并点名出现
  的负年龄与父子倒置。0.1.0 遇到同样的校正表，会把结果连着警告一起返回。本教程
  此前那一版占位示例——只有密度的 `Elephas`/`Loxodonta` 行，配上 `70`–`105` 的
  `Dasypus`/`Choloepus` 边界——正是这样一种不可满足的表；随包的
  `data/examples/example_calibrations.tsv` 也仍然是，所以它被写成一个"预期会失败"
  的格式示意，见 `data/examples/README-zh.md`。想拿到开箱即用的表，请用
  `data/examples/example_calibrations_feasible.tsv`——同样的四条硬边界，只是没有
  密度那一行。

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
