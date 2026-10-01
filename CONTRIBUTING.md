# Contributing to OpenRelTime

Thanks for your interest!  This document covers what a patch needs to get
accepted.  中文版见 [CONTRIBUTING-zh.md](CONTRIBUTING-zh.md)。Issues 与 PR
使用英文或中文均可。

## Development setup

```bash
git clone https://github.com/ZengZichao/OpenRelTime
cd OpenRelTime
python -m pip install --upgrade pip
python -m pip install --editable ".[dev,plot]"
# dev = pytest, pytest-cov, ruff, mypy;  plot = matplotlib
```

Optional externals: the `blb` and `megacc` end-to-end cases shell out to
**IQ-TREE** and **MEGA-CC**.  Without them installed, those two modules
skip; everything else runs.

## Checks every PR must pass

```bash
python -m ruff check . --select F,E9     # syntax / unused names
python -m pytest                         # unit suite (tests/)
python -m pytest e2e_tests/cases         # end-to-end suite
python -m build                          # only when packaging files change
```

CI runs the same checks on Ubuntu (Python 3.10–3.13) and macOS (3.13), plus
a `minimum-versions` job that pins the oldest declared dependency floors
(`.github/constraints-minimum.txt`).  `mypy` currently runs **advisory** in
CI — run it locally with
`python -m mypy openreltime --python-version 3.12` and avoid introducing new
findings.

## Conventions

- Commit messages follow Conventional Commits
  (`feat:`, `fix:`, `docs:`, `chore:`, ...).
- User-facing documentation is bilingual: if you change `docs/*.md` or a
  README, update the `-zh` twin in the same PR.
- Golden data under `data/golden/` and `e2e_tests/data/golden/` is
  regenerated with the scripts in `reproduce/` — never edited by hand.
  If a change legitimately alters outputs, regenerate the golden files and
  explain why in the PR description.

## Submitting

Open a pull request against `main` using the PR template.  The `all-green`
check must pass before merge.
