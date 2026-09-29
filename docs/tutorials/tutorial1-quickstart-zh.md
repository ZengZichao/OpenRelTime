# 教程 1：五分钟上手 OpenRelTime

**语言：中文** · [English](tutorial1-quickstart.md)

```python
import openreltime as ort

# 1. 读入带枝长的树；以外群置根并剔除外群
tree = ort.read_tree(
    "data/examples/example.nwk",
    outgroup=[
        "Ornithorhynchus_anatinus",
        "Zaglossus_bruijni",
        "Tachyglossus_aculeatus",
    ],
)

# 2. 谱系相对速率与节点相对时间
rates = ort.rrf_rates(tree)
times = ort.rrf_rates_times(tree)          # 一次算齐两者

# 3. 写出 CSV 表、时间树与带速率注释的 NEXUS
times.write("tutorial1", with_rate=True, nexus=True)

# 4. 分子钟是否存在速率自相关？
corr = ort.corrtest(tree)
print(corr.score, corr.p_band)

# 5. 由相对时间估计出生-死亡物种分化先验
bd = ort.ddbd(tree, anchor_time=1.85)
print(bd.birth_rate, bd.death_rate, bd.sampling_frac)
```

在 shell 中完成同样的分析：

```bash
openreltime rates-times -i example.nwk --outgroup "Ornithorhynchus_anatinus,Zaglossus_bruijni,Tachyglossus_aculeatus" -o tutorial1 --plot tutorial1.png
openreltime corrtest -i example.nwk --outgroup "Ornithorhynchus_anatinus,Zaglossus_bruijni,Tachyglossus_aculeatus" -o tutorial1
openreltime ddbd -i example.nwk --outgroup "Ornithorhynchus_anatinus,Zaglossus_bruijni,Tachyglossus_aculeatus" --anchor-time 1.85 -o tutorial1
```

写出的文件如下表。三者都共用 `tutorial1` 前缀，三条命令会把各自的参数、警告与随
机种子合并写入同一个 `tutorial1_report.json`。

| 文件 | 来自 |
|---|---|
| `tutorial1_rates_times.csv` | `rates-times`（速率与时间两列；加 `--r3f-compat` 会改名为 `_RRF_table.csv`） |
| `tutorial1_timetree.nexus` | `rates-times`（带 `[&rate=...]` 注释的相对时间树，FigTree 可直接打开） |
| `tutorial1.png` | `rates-times --plot`（按速率着色的时间树图像） |
| `tutorial1_corrtest.txt` | `corrtest`（CorrScore 与 P 值档位） |
| `tutorial1_ddbd.txt` | `ddbd`（出生、死亡与抽样比例） |
| `tutorial1_report.json` | 每一次运行（参数、警告与随机种子，用于复现） |

注意：`rates-times` 写出 `_rates_times.csv` 表与 `_timetree.nexus` 树，**不**额外
写出单独的 `_rates.csv` 和 `_timetree.nwk`。需要独立的速率表，请用 `rates` 子命
令（输出 `<prefix>_rates.csv` 与 `<prefix>_rates.nwk`）；需要普通的 Newick 时间
树，请用 `times` 子命令（输出 `<prefix>_times.csv` 与
`<prefix>_timetree.nwk`）。

下一步：[教程 2（端到端定年）](tutorial2-end-to-end-zh.md)。
