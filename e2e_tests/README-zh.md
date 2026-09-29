# OpenRelTime 端到端测试套件

**语言：中文** · [English](README.md)

本目录包含 OpenRelTime 的**全流程真实数据端到端测试**。这些测试不模拟内部函数，而是通过命令行（CLI）与 Python API 真实运行 `openreltime` 的各个分析模块，并与项目自带的 golden 结果进行比对。

测试文件随软件包一起发布，用户可在安装后独立运行，以验证本地环境能够完整复现示例分析流程。

## 目录结构

```
e2e_tests/
├── README.md                 # English readme
├── README-zh.md              # 本说明
├── conftest.py               # 共享 fixtures、路径、CLI 运行器
├── data/
│   ├── example.nwk           # 示例系统树（Newick）
│   ├── example_taxa.tsv      # 示例分组表
│   ├── example_calibrations.tsv  # 示例校正点
│   ├── example.fasta         # 较大的随机生成序列对齐（可选）
│   ├── small.fasta           # 用于 BLB 冒烟测试的极小随机生成对齐
│   ├── golden/               # golden 结果（复制自 data/golden/r3f/）
│   └── bad/                  # 错误输入用例的临时目录（当前为空，测试按需生成）
├── cases/
│   ├── test_01_rrf_core.py
│   ├── test_02_calibration_ci.py
│   ├── test_03_corrtest_ddbd.py
│   ├── test_04_tree2table_monophyly.py
│   ├── test_05_cli_api_parity.py
│   ├── test_06_input_validation.py
│   ├── test_07_blb_external.py
│   └── test_08_megacc_external.py
└── outputs/                  # 测试输出目录（被 .gitignore 忽略）
```

## 环境要求

1. 激活 `openreltime` conda 环境：
   ```bash
   conda activate openreltime
   ```
2. 以 editable 模式安装当前包及全部可选依赖：
   ```bash
   pip install --editable ".[dev,plot]"
   ```
3. （可选）若希望运行外部工具测试，需安装对应二进制程序：
   - **IQ-TREE**：用于 BLB 测试。可用 micromamba 安装到当前环境：
     ```bash
     micromamba install -p /path/to/openreltime/env -c conda-forge -c bioconda iqtree
     ```
   - **MEGA-CC `megacc`**：用于 megacc 桥接测试。该工具未在 conda 主频道提供，需手动安装并将路径写入环境变量 `MEGACC`。

## 测试数据说明

| 文件 | 来源 | 用途 |
|------|------|------|
| `data/example.nwk` | `data/examples/example.nwk` | RRF、校正、CI、CorrTest、ddBD、monophyly 等所有树分析 |
| `data/example_taxa.tsv` | `data/examples/example_taxa.tsv` | monophyly / tree2table 分组测试 |
| `data/example_calibrations.tsv` | `data/examples/example_calibrations.tsv` | calibrate、CI 与 megacc 测试 |
| `data/small.fasta` | 随机生成（10 条序列 × 300 bp） | BLB / IQ-TREE 外部流程冒烟测试 |
| `data/example.fasta` | 随机生成（274 条） | 保留作为可选大样本对齐 |
| `data/golden/*` | `data/golden/r3f/` | 与 CLI/API 输出进行回归比对 |

所有测试数据均位于 `e2e_tests/data/` 内，测试过程中不会在仓库根目录之外创建或修改文件。

## 运行方式

```bash
# 完整套件（推荐）
python -m pytest e2e_tests/ -v

# 仅运行 RRF 核心测试
python -m pytest e2e_tests/cases/test_01_rrf_core.py -v

# 生成 HTML 报告
python -m pytest e2e_tests/ -v --html=e2e_tests/outputs/report.html
```

> 说明：`pyproject.toml` 中的默认 `testpaths` 只包含 `tests/`，因此运行 e2e
> 套件时需要显式指定 `e2e_tests/` 路径。

## 测试模块与覆盖范围

| 模块 | 用例数 | 覆盖功能 |
|------|--------|----------|
| `test_01_rrf_core.py` | 7 | `rates` / `times` CLI、golden 回归、输出文件、归一化、算术平均、guard 选项、NEXUS 往返 |
| `test_02_calibration_ci.py` | 9 | `calibrate` bounds/effective、CI 默认/branch-var/n-sites、API 与 CLI 一致性、冲突校正、纯密度约束、pinned node 标记 |
| `test_03_corrtest_ddbd.py` | 7 | CorrTest golden 回归、全精度 API、重采样；ddBD golden、anchor、plot、API/CLI 一致性 |
| `test_04_tree2table_monophyly.py` | 5 | `tree2table` 支长与节点年龄；monophyly 自动/表格/临时分组 |
| `test_05_cli_api_parity.py` | 3 | API 写出的结果与 CLI 一致、`_report.json` 合并、参数记录 |
| `test_06_input_validation.py` | 8 | 负支长、重复 tip、多岔、外类群不完整、空校正行、密度-only 与方法冲突、缺失输入文件 |
| `test_07_blb_external.py` | 2 | BLB 全流程（IQ-TREE）与自定义 IQ-TREE 参数传递 |
| `test_08_megacc_external.py` | 2 | MEGA-CC 桥接：含校正点与不含校正点 |
| **合计** | **43** | 8 个模块，仅覆盖引擎 CLI 与 API |

## 跳过条件

- `test_07_blb_external.py`：当 `IQTREE` 环境变量未设置且系统 PATH 中找不到 `iqtree2` / `iqtree` 时自动跳过。
- `test_08_megacc_external.py`：当 `MEGACC` 环境变量未设置且 PATH 中找不到 `megacc` 时自动跳过。

## 输出目录

测试生成的所有文件均写入 `e2e_tests/outputs/`，该目录已被加入 `.gitignore`，不会进入版本控制。每次测试的单个用例输出位于 `e2e_tests/outputs/<pytest_tmp_name>/` 下。

## 注意事项

- 测试使用 `subprocess.run` 调用 `python -m openreltime.cli`，确保验证的是用户实际使用的 CLI 入口。
- 部分用例涉及随机数，已通过 `--seed` 固定，保证结果可复现。
- 本套件只覆盖引擎。桌面图形界面作为独立产品 `OpenRelTime-Studio` 发行，测试
  套件也在该项目内；本仓库的测试既不导入也不启动它。
