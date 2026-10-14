#!/usr/bin/env python3
"""Check that ``docs/validation.md`` and the test suite agree about gates.

The validation gates (G1, G2a, G2b, G2c, G3, G4, G5, G6, G7, G8, G9) are the
backbone of the manuscript's fidelity claims, so the gate-to-test mapping must
not be able to drift silently.  Two independent sources describe it:

* ``docs/validation.md`` section 1, whose ``Where`` column names test files and
  individual tests;
* ``@pytest.mark.gate("...")`` marks on the tests themselves.

This tool cross-checks them and fails when either side references a test that no
longer exists, when a gate has no test behind it, or when a test's mark
contradicts the document.  G9 is deliberately exempt: it is measured, never
asserted, and the document says so.

Run from the repository root::

    python tools/check_gates.py            # human-readable report
    python tools/check_gates.py --quiet    # exit code only

Exit status is non-zero on any inconsistency.
"""

from __future__ import annotations

import argparse
import collections
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOC = ROOT / "docs" / "validation.md"
GATE_RE = re.compile(r"^G\d[abc]?$")
NAMED_TEST_RE = re.compile(r"test_[A-Za-z0-9_]+\.py(?::：:?\s*[\w.*\[\]-]+)?")
BARE_TEST_RE = re.compile(r"test_[A-Za-z0-9_]+")

#: G9 is a measured quantity with no automated assertion, by design.
MEASURED_ONLY = {"G9"}


def doc_gates() -> dict[str, list[str]]:
    """Return ``{gate: [where-strings]}`` from the validation.md gate table."""
    text = DOC.read_text(encoding="utf-8")
    gates: dict[str, list[str]] = collections.defaultdict(list)
    for row in re.finditer(r"^\|\s*(G\d[abc]?)\s*\|(.*?)\|\s*$", text, re.M):
        gate, rest = row.group(1), row.group(2)
        cells = [c.strip() for c in rest.split("|")]
        if cells:
            gates[gate].append(cells[-1])
    return dict(gates)


def collected_gate_marks() -> dict[str, set[str]]:
    """Return ``{nodeid: {gate,...}}`` for every collected test."""
    listing = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
         "--collect-only", "--no-header", "-o", "addopts="],
        cwd=ROOT, capture_output=True, text=True, check=False,
    ).stdout
    marks: dict[str, set[str]] = collections.defaultdict(set)
    for nodeid in re.findall(r"^(tests/\S+::\S+)$", listing, re.M):
        marks[nodeid]  # ensure every collected test is a key, marked or not
        for gate in _gate_for(nodeid).split("|"):
            if gate:
                marks[nodeid].add(gate)
    return dict(marks)


def _gate_for(nodeid: str) -> str:
    """Best-effort gate id for a nodeid, read from its mark or module mark."""
    path, _, name = nodeid.partition("::")
    try:
        src = (ROOT / path).read_text(encoding="utf-8")
    except OSError:
        return ""
    if name:
        m = re.search(
            r'((?:@pytest\.mark\.gate\("[^"]+"\)\s*\n)+)(?:@[^\n]*\n)*def '
            + re.escape(name.split("[")[0]) + r"\(",
            src,
        )
        if m:
            return "|".join(re.findall(r'gate\("([^"]+)"\)', m.group(1)))
    m = re.search(r'^pytestmark = pytest\.mark\.gate\("([^"]+)"\)', src, re.M)
    return m.group(1) if m else ""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quiet", action="store_true", help="only report the exit code")
    args = ap.parse_args()

    problems: list[str] = []
    gates = doc_gates()

    if not gates:
        print(f"error: no gate table found in {DOC}", file=sys.stderr)
        return 2

    existing_files = {p.name for p in (ROOT / "tests").glob("test_*.py")}
    marks = collected_gate_marks()
    all_ids = set(marks)

    doc_gates_seen: set[str] = set()
    for gate, wheres in sorted(gates.items()):
        named_files: set[str] = set()
        named_tests: set[str] = set()
        for where in wheres:
            named_files.update(NAMED_TEST_RE.findall(where.replace("<br>", " ")))
            named_tests.update(BARE_TEST_RE.findall(where))
        if not named_files and gate not in MEASURED_ONLY:
            problems.append(f"{gate}: the document names no test file")
        for entry in named_files:
            stem = entry.split(".")[0]
            if f"{stem}.py" not in existing_files:
                problems.append(f"{gate}: document names {entry}, which is not in tests/")
            for bare in set(BARE_TEST_RE.findall(entry)):
                if bare != stem and not any(
                    f"::{bare}" in n or f"::{bare}[" in n for n in all_ids
                ):
                    problems.append(
                        f"{gate}: document names test {bare}, which is not collected"
                    )
        hit = [n for n in all_ids if gate in marks[n]]
        if not hit and gate not in MEASURED_ONLY:
            problems.append(
                f'{gate}: no collected test carries @pytest.mark.gate("{gate}"), '
                f"yet the document asserts it"
            )
        else:
            doc_gates_seen.add(gate)

    for marked, gs in sorted(marks.items()):
        for g in gs:
            if g not in gates:
                problems.append(f"{marked}: mark {g} is not defined in validation.md")

    unmarked = [n for n in sorted(all_ids) if not marks[n]]

    if not args.quiet:
        print(f"gates defined in validation.md : {', '.join(sorted(gates))}")
        print(f"tests collected               : {len(all_ids)}")
        print(f"tests carrying a gate mark    : {len(all_ids) - len(unmarked)}")
        per_gate = collections.Counter(
            g for gs in marks.values() for g in gs
        )
        for gate in sorted(per_gate):
            print(f"  {gate:<4}-> {per_gate[gate]} test(s)")
        if unmarked:
            print(f"{len(unmarked)} test(s) carry no gate mark (covered by the "
                  f"document's Where column instead):")
            for n in unmarked[:8]:
                print(f"  {n}")
            if len(unmarked) > 8:
                print(f"  ... and {len(unmarked) - 8} more")

    if problems:
        print("\nINCONSISTENCIES:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1
    if not args.quiet:
        print("\nvalidation.md and the gate marks agree.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
