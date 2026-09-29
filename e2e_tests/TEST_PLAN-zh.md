# OpenRelTime 端到端测试方案

**语言：中文** · [English](TEST_PLAN.md)

## 1. 测试目标

本端到端（End-to-End, E2E）测试方案旨在使用**真实示例数据**对 OpenRelTime 的完整分析流程进行自动化验证，覆盖从命令行（CLI）到 Python API、从核心算法到外部工具桥接、从输入校验到输出可复现性的全部主要功能。所有测试文件位于仓库 `e2e_tests/` 目录下，随源码包一起发布，确保用户在任何环境中都能复现验证。

## 2. 测试范围

| 功能域 | 覆盖内容 |
|--------|----------|
| RRF 核心 | `rates`、`times` CLI；golden 回归；输出文件（CSV、NEXUS、Newick）；归一化；算术平均；guard 选项；NEXUS 往返 |
| 校正与置信区间 | `calibrate`（bounds/effective）；`ci`（默认、branch-var、n-sites）；API 与 CLI 一致性；冲突校正；纯密度约束；pinned node 标记 |
| 相关性 / 出生死亡 | `corrtest` golden 回归、全精度 API、重采样；`ddbd` golden、anchor、plot、API/CLI 一致性 |
| 树表与单系性 | `tree2table` 支长与节点年龄；`monophyly` 自动/表格/临时分组 |
| CLI/API 等价性 | API 写出结果与 CLI 输出一致、`_report.json` 合并、参数记录 |
| 输入校验 | 负支长、重复 tip、多岔、外类群不完整、空校正行、密度-only 与方法冲突、缺失输入文件 |
| 外部工具 | BLB（IQ-TREE）；MEGA-CC（megacc）桥接 |

## 3. 测试数据

所有数据均位于 `e2e_tests/data/`，不超出工作路径。

| 文件 | 来源/生成方式 | 用途 |
|------|--------------|------|
| `example.nwk` | `data/examples/example.nwk` | 主系统树 |
| `example_taxa.tsv` | `data/examples/example_taxa.tsv` | 分组表 |
| `example_calibrations.tsv` | `data/examples/example_calibrations.tsv` | 校正点 |
| `small.fasta` | 随机生成（10 条 × 300 bp） | BLB / IQ-TREE 冒烟测试 |
| `example.fasta` | 随机生成（274 条） | 保留的可选大样本对齐 |
| `golden/` | `data/golden/r3f/` | 回归比对的基准结果 |

## 4. 测试环境

- **操作系统**：Linux x86_64
- **Python**：≥3.10
- **核心依赖**：numpy、scipy、pandas、click
- **可选依赖**：matplotlib、pytest
- **外部工具**：IQ-TREE（BLB 测试）、MEGA-CC `megacc`（megacc 桥接测试）
- **安装命令**：
  ```bash
  pip install --editable ".[dev,plot]"
  micromamba install -p /path/to/openreltime/env -c conda-forge -c bioconda iqtree
  ```

## 5. 测试执行策略

- 使用 `pytest` 作为测试框架，测试文件统一命名为 `test_*.py`。
- 通过 `subprocess.run` 调用 `python -m openreltime.cli`，验证真实 CLI 行为。
- 外部工具测试通过 `skipif` 自动检测可执行文件；缺失时跳过并记录。
- 每个测试使用独立的 `tmp_out` 目录，避免相互污染。

## 6. 用例设计与验证标准

### 6.1 RRF 核心（`test_01_rrf_core.py`）

| 用例 | 输入 | 验证标准 |
|------|------|----------|
| `test_rates_cli_regression` | `example.nwk` + 外类群 | CLI 生成 `_rates.csv`、`_rates.nwk`，关键节点 rate 与 golden 一致 |
| `test_times_cli_r3f_compat` | `example.nwk` + 外类群 | 生成 `_times.csv`、`_timetree.nwk`，根节点时间与 golden 一致 |
| `test_rates_times_cli_outputs` | 同上 | 同时运行 `rates-times`，`_report.json` 包含 rates/times 两个分析块 |
| `test_times_normalize` | `--normalize` | 节点时间整体缩放，结构不变 |
| `test_arithmetic_mean_cli` | `--mean arithmetic` | 使用算术平均成功运行 |
| `test_guard_options` | `--no-guard` / 默认 | guard 关闭时根节点时间不同于默认 |
| `test_nexus_roundtrip` | NEXUS 输出再读入 | 读写 NEXUS 不丢失节点数 |

### 6.2 校正与 CI（`test_02_calibration_ci.py`）

| 用例 | 输入 | 验证标准 |
|------|------|----------|
| `test_calibrate_bounds_cli` | 校正点 + `--method bounds` | 输出 `_calibrated_*.csv`、time factor > 0 |
| `test_calibrate_effective_cli` | 校正点 + `--method effective` | effective 方法成功运行 |
| `test_ci_default` | 已校正结果 | 输出 CI 文件、节点年龄上下界 |
| `test_ci_with_n_sites` | `--n-sites` | CI 宽度随 n-sites 变化 |
| `test_ci_with_branch_var` | `--branch-var` | 使用支长方差运行 |
| `test_calibration_api_matches_cli` | API 与 CLI 同输入 | 两者 time factor 一致 |
| `test_conflicting_calibration_rejected` | 边界冲突的校正点 | CLI 返回非 0 并提示冲突 |
| `test_density_only_requires_effective` | 密度-only 校正 + effective | 成功运行 |
| `test_pinned_node_ci_flagged` | pinned 校正点 | CI 表中对应节点标记为 pinned |

### 6.3 CorrTest / ddBD（`test_03_corrtest_ddbd.py`）

| 用例 | 输入 | 验证标准 |
|------|------|----------|
| `test_corrtest_cli_matches_golden` | `example.nwk` | CorrScore 与 golden 一致 |
| `test_corrtest_api_full_precision` | 同上 | API 全精度输出与 golden 一致 |
| `test_corrtest_resample` | `--sister-resample 50` | 重采样后结果与 `sr50` golden 一致 |
| `test_ddbd_cli_matches_golden` | 同上 | birth/death rate 与 golden 一致 |
| `test_ddbd_anchor_matches_golden` | `--anchor-node` / `--anchor-time` | anchor 结果与 golden 一致 |
| `test_ddbd_plot` | 同上 | 生成 `_ddbd.png` |
| `test_ddbd_api_matches_cli` | API 与 CLI 同输入 | 两者数值一致 |

### 6.4 tree2table / monophyly（`test_04_tree2table_monophyly.py`）

| 用例 | 输入 | 验证标准 |
|------|------|----------|
| `test_tree2table_branch_lengths` | `example.nwk` | CSV 包含原始支长 |
| `test_tree2table_node_ages` | RRF 时间树 | CSV 节点年龄非空 |
| `test_monophyly_auto_groups` | `--taxa-table` | 自动分组结果符合预期 |
| `test_monophyly_taxon_table_groups` | 同上 | 表格中各组单系性判断正确 |
| `test_monophyly_ad_hoc_tips` | `--tips` | 临时 tip 集合判断正确 |

### 6.5 CLI/API 等价性（`test_05_cli_api_parity.py`）

| 用例 | 输入 | 验证标准 |
|------|------|----------|
| `test_api_write_matches_cli_rates_times` | API 生成文件与 CLI 对比 | DataFrame 差值为 0 |
| `test_report_json_merges_blocks` | 多阶段分析 | `_report.json` 合并所有阶段参数 |
| `test_cli_parameters_recorded` | 带参数运行 | report 中参数与命令行一致 |

### 6.6 输入校验（`test_06_input_validation.py`）

| 用例 | 输入 | 验证标准 |
|------|------|----------|
| `test_negative_branch_length_rejected` | 负支长树 | 返回非 0，stderr 含 "negative" |
| `test_duplicate_tip_rejected` | 重复 tip 树 | 返回非 0，stderr 含 "duplicate" |
| `test_polytomy_rejected_by_default` | 多岔树 | 默认报错，`--resolve random` 可通过 |
| `test_too_few_ingroup_tips_rejected` | 外类群包含内类群 | 返回非 0 |
| `test_incomplete_outgroup_rejected` | 不完整外类群 | 默认报错，`--outgroup-check warn` 可通过 |
| `test_empty_calibration_row_rejected` | 空校正行 | 返回非 0，stderr 含 "bound" |
| `test_density_only_with_bounds_method_rejected` | 密度-only + bounds | 返回非 0，stderr 含 "effective" |
| `test_missing_input_file` | 不存在的输入 | 返回非 0，stderr 含 "not found" |

### 6.7 外部工具（`test_07_blb_external.py`、`test_08_megacc_external.py`）

| 用例 | 输入 | 验证标准 |
|------|------|----------|
| `test_blb_full_pipeline` | `small.fasta` + IQ-TREE | 生成 `blb_summary.csv` 与 `blb_replicates.csv` |
| `test_blb_iqtree_args` | `--iqtree-arg` | 自定义参数不报错 |
| `test_megacc_with_calibrations` | `example.nwk` + 校正点 + megacc | 生成 `_megacc_times.csv` |
| `test_megacc_without_calibrations` | `example.nwk` + megacc | 生成 `_megacc_times.csv` |

## 7. 跳过条件

- **IQ-TREE 测试**：当 `IQTREE` 环境变量未设置且 PATH 中不存在 `iqtree2` / `iqtree` 时跳过。
- **MEGA-CC 测试**：当 `MEGACC` 环境变量未设置且 PATH 中不存在 `megacc` 时跳过。

## 8. 风险与假设

1. 外部工具版本差异可能导致数值微小漂移；golden 测试使用固定容差。
2. IQ-TREE 的 `-nt AUTO` 在部分机器上极慢，因此 BLB 冒烟测试固定 `-nt 1 -m GTR+G` 并使用最小对齐。
3. MEGA-CC 未在 conda 主频道提供，视为可选外部依赖。
4. 桌面图形界面不在本方案范围内：它作为独立产品 `OpenRelTime-Studio` 发行，
   测试也在该项目内完成。本方案只检验引擎的 CLI 与 Python API。

## 9. 回归与维护

- 新增主要功能时，应在 `e2e_tests/cases/` 下新增对应 `test_*.py`。
- 修改 CLI 参数或输出格式时，必须同步更新 e2e 测试与 golden 数据。
- 每次发布前执行：
  ```bash
  python -m pytest e2e_tests/ -v
  ```
