#!/usr/bin/env python
"""Print the collected test count of every suite in this repository.

The manuscript and ``docs/validation.md`` quote per-suite test counts, and those
numbers drift silently whenever a test is added.  This script makes them
machine-checkable: it runs ``pytest --collect-only`` for each suite and prints
the collected count, per file, optionally as JSON for CI or doc generation.

Suites (they are *not* all collected by a bare ``python -m pytest``)::

    tests/       the engine unit suite; this is ``[tool.pytest.ini_options]``
                 ``testpaths``, so it is what CI's "Run test suite" step runs
    e2e_tests/   end-to-end CLI suite; CI runs it as a separate step
                 (``python -m pytest e2e_tests/cases``), and it is excluded from
                 the default ``testpaths`` on purpose -- adding it there would
                 change what the default invocation (and CI) means

Usage::

    python tools/count_tests.py            # human-readable table
    python tools/count_tests.py --json     # machine-readable
    python tools/count_tests.py -k G7      # only tests carrying that gate mark

Exit code is 0 when every suite collected, 1 when any suite failed to collect.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]

#: name -> (path passed to pytest, extra argv).  ``-c /dev/null`` neutralises
#: ``testpaths`` so a named path is collected even though the default
#: invocation would only ever look in ``tests/``.
SUITES: dict[str, tuple[str, list[str]]] = {
    "tests": ("tests", []),
    "e2e_tests": ("e2e_tests", ["-c", "/dev/null"]),
}

_COLLECTED = re.compile(r"(\d+)\s+tests?\s+collected")


def _collect(path: str, extra: list[str], keyword: str | None) -> tuple[int, list[str]]:
    """Return (collected count, per-file "path::count" lines) for one suite."""
    argv = [
        sys.executable,
        "-m",
        "pytest",
        "--collect-only",
        "-q",
        "-p",
        "no:cacheprovider",
        "-o",
        "addopts=",
        *extra,
    ]
    if keyword:
        argv += ["-k", keyword]
    argv.append(path)
    proc = subprocess.run(
        argv, cwd=str(ROOT), capture_output=True, text=True, check=False
    )
    text = proc.stdout + proc.stderr
    match = _COLLECTED.search(text)
    if match is None:
        if "no tests ran" in text or proc.returncode in (5,):
            return 0, []
        raise RuntimeError(
            f"could not read a collected-test count for {path!r}\n{text[-800:]}"
        )
    lines = [
        line
        for line in proc.stdout.splitlines()
        if line.strip().endswith(".py") or "::" in line
    ]
    per_file: dict[str, int] = {}
    for line in lines:
        file_part = line.split("::")[0].strip()
        if file_part.endswith(".py"):
            per_file[file_part] = per_file.get(file_part, 0) + 1
    return int(match.group(1)), [
        f"{name}::{count}" for name, count in sorted(per_file.items())
    ]


def main(argv: list[str] | None = None) -> int:
    """Count every suite and print it; return the process exit code."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--json", action="store_true", help="emit JSON")
    parser.add_argument(
        "-k",
        "--keyword",
        default=None,
        help="pytest -k expression, e.g. a gate id, to subset the count",
    )
    args = parser.parse_args(argv)

    result: dict[str, Any] = {}
    failed = False
    for name, (path, extra) in SUITES.items():
        try:
            count, files = _collect(path, extra, args.keyword)
        except RuntimeError as exc:
            print(f"{name}: ERROR {exc}", file=sys.stderr)
            failed = True
            count, files = -1, []
        result[name] = {"path": path, "collected": count, "files": files}
    result["total"] = sum(
        entry["collected"] for entry in result.values() if isinstance(entry, dict)
    )

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        width = max(len(name) for name in SUITES) + 2
        for name, entry in result.items():
            if name == "total":
                continue
            print(f"{name:<{width}s}{entry['collected']:>6d}  ({entry['path']})")
        print(f"{'-' * (width + 6)}")
        print(f"{'total':<{width}s}{result['total']:>6d}")
        print(
            "\nNOTE: a bare `python -m pytest` collects only tests/ (testpaths in "
            "pyproject.toml);\ne2e_tests/ is a separate CI step. See "
            "docs/validation.md."
        )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
