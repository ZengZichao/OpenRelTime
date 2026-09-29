# 教程 3：用 ddBD 生成 MCMCTree 树先验，以及 MEGA-CC 桥接

**语言：中文** · [English](tutorial3-priors-megacc.md)

## 供贝叶斯定年使用的出生-死亡先验

```python
import openreltime as ort

tree = ort.read_tree("mam.treefile", outgroup=["Out1", "Out2"])

# 指定要报告的抽样比例（出生率/死亡率仍来自自由 rho 的拟合）
bd = ort.ddbd(tree, sampling_frac=0.5, anchor_time=1.85)

# 全部待估（默认）
bd_full = ort.ddbd(tree, anchor_time=1.85, measure="KL")

print(bd.birth_rate, bd.death_rate, bd.sampling_frac)
bd.write("mam_ddbd")
```

估计得到的出生率、死亡率与抽样比例共同参数化一个出生-死亡节点年龄密度，可直接
用作 MCMCTree 与 BEAST 的树先验。请保持 `anchor_time`（以及由此决定的最大节点年
龄）≤ 10 个时间单位。

`sampling_frac` 的语义与 R3F 一致：该参数设定*报告*的抽样比例。似然优化始终把
rho 留作自由参数，因此 `sampling_frac` 并不约束拟合。自由拟合得到的 rho 仍可在
报告的 `sampling_frac_fitted` 中查看。

## MEGA-CC 桥接（可选 extra）

供从 MEGA 流程迁移或做第三方对照使用：

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

`run_megacc` 生成一个最小的 `reltimeFromBranchLengths.mao`
（`ppRelTimeBLens = true` 加上可选覆盖项），探测 megacc 的版本横幅，并在启动进程
之前校验每一个参数（不经 shell）。程序把 megacc 的输出解析回与 R3F 同构的节点时
间表，`source` 列记录每一行的来源，因此 MEGA 结果可与 `openreltime.rrf_times` 的
输出逐节点比较。

对应命令行：

```bash
openreltime megacc -i mam.treefile -c calibrations.tsv --megacc /opt/megacc \
    --outgroup "Out1,Out2" --workdir megacc_run -o mam_mega
# 结果写入 mam_mega_megacc_times.csv
```

桥接有两个限制，详见 `docs/validation-zh.md`。第一，它需要真实的 `megacc` 可执行
文件，CI 既不安装也不运行该程序，因此 MEGA-CC 的一致性依据已发表的 R3F/MEGA 结果与
随包黄金文件确立，而不是依据此处新做的头对头运行。第二，OpenRelTime
的解析 RRF 层缺少 MEGA 的似然比精修，逐节点差异属预期之内。把它当作对照工具，
不要当作等价性声明。

上一步：[教程 2](tutorial2-end-to-end-zh.md)。
