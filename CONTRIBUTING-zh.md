# 参与贡献 OpenRelTime

感谢关注！本文档说明一个补丁被接受所需满足的最低要求。
English version: [CONTRIBUTING.md](CONTRIBUTING.md)。Issue 与 PR 使用中文
或英文均可。

## 开发环境

```bash
git clone https://github.com/ZengZichao/OpenRelTime
cd OpenRelTime
python -m pip install --upgrade pip
python -m pip install --editable ".[dev,plot]"
# dev = pytest、pytest-cov、ruff、mypy；plot = matplotlib
```

可选外部依赖：`blb` 与 `megacc` 两个端到端用例会调用 **IQ-TREE** 和
**MEGA-CC**。未安装时这两个模块会跳过，其余用例正常运行。

## 每个 PR 必须通过的检查

```bash
python -m ruff check . --select F,E9     # 语法 / 未使用名称
python -m pytest                         # 单元测试（tests/）
python -m pytest e2e_tests/cases         # 端到端测试
python -m build                          # 仅在改动打包相关文件时
```

CI 在 Ubuntu（Python 3.10–3.13）与 macOS（3.13）上运行同样的检查，另有
`minimum-versions` 作业按 `.github/constraints-minimum.txt` 钉住依赖的
最低声明版本。`mypy` 目前在 CI 中为**咨询模式**（不阻塞），本地可运行
`python -m mypy openreltime --python-version 3.12`，请避免引入新的类型
问题。

## 约定

- 提交信息遵循 Conventional Commits
  （`feat:`、`fix:`、`docs:`、`chore:` 等）。
- 面向用户的文档为中英双语：改动 `docs/*.md` 或任一 README 时，请在
  同一个 PR 中同步更新对应的 `-zh` 版本。
- `data/golden/` 与 `e2e_tests/data/golden/` 下的 golden 数据由
  `reproduce/` 中的脚本再生成，**不要手改**。若改动确实会改变输出，
  请再生成 golden 文件并在 PR 描述中说明原因。

## 提交方式

按 PR 模板向 `main` 发起 pull request。合并前 `all-green` 检查必须通过。
