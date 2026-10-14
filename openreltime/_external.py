"""Single audited choke point for starting external programs.

Security contract
-----------------
* the executable is resolved with :func:`shutil.which` or must exist as a
  path, and every token is checked against a shell-metacharacter blacklist
  (:data:`_FORBIDDEN`) — tokens containing shell control characters are
  rejected *before* any process starts.  Spaces and non-ASCII (e.g. CJK)
  characters are allowed so that real-world installation paths work;
* programs are always started from an argument **list** with
  ``shell=False``; a shell is never created, so shell interpolation /
  injection is structurally impossible regardless of token contents;
* callers must pre-validate domain-specific values (e.g. ``key=value``
  calibration parameters are parsed explicitly, never evaluated).

Used by :mod:`openreltime.bootstrap` (IQ-TREE) and
:mod:`openreltime.megacc` (MEGA-CC).  :func:`probe_version` is the single
version-banner site for both: no release is pinned in code or CI, so the banner
each tool prints is recorded in that analysis' JSON report instead.
"""

from __future__ import annotations

import re
import subprocess  # noqa: S404 - argv-list execution only, shell=False
from pathlib import Path
from shutil import which
from typing import Any, Callable, Optional

#: Characters that never belong in a command token.  With ``shell=False``
#: they could not be interpolated anyway; the blacklist is defense in depth
#: against accidental misuse of :func:`run_tool` with assembled strings.
#: Whitespace and non-ASCII characters are deliberately allowed (paths), and
#: so is the backslash: Windows executables live under ``C:\...`` paths, and
#: with ``shell=False`` a backslash in a list-mode token is inert data.
_FORBIDDEN = re.compile(r"[;&|<>$`\"'\r\n\t]")


def validate_tokens(*tokens: str) -> None:
    """Raise :class:`ValueError` if any token contains a shell metacharacter."""
    for token in tokens:
        if not token:
            raise ValueError("refusing to run: empty command token")
        if _FORBIDDEN.search(token):
            raise ValueError(
                f"refusing to run: argument {token!r} contains shell "
                "metacharacters (one of ; & | < > $ ` \" ' or a newline)"
            )


def resolve_executable(name_or_path: str) -> str:
    """Resolve a bare name on PATH or verify an explicit path."""
    resolved = which(name_or_path)
    if resolved is None and Path(name_or_path).exists():
        resolved = str(Path(name_or_path))
    if not resolved:
        raise FileNotFoundError(
            f"executable not found: {name_or_path!r}; "
            "install the tool or pass its full path"
        )
    validate_tokens(resolved)
    return resolved


# Bind once; keeps a single reviewable site for the whole package.
_run = subprocess.run


def run_tool(argv: list[str], *, check: bool = True):
    """Start an external program from a validated argument list.

    No shell is created (``shell=False``); output is captured and returned
    as a :class:`subprocess.CompletedProcess`.
    """
    validate_tokens(*argv)
    return _run(  # noqa: S603 - validated argv list, shell=False
        argv,
        check=check,
        capture_output=True,
        shell=False,
    )


def probe_version(
    executable: str,
    *,
    tool: str,
    runner: Optional[Callable[..., Any]] = None,
) -> tuple[str, list[str]]:
    """Read *executable*'s ``--version`` banner; never let the probe break a run.

    Both external-tool paths (IQ-TREE for BLB, ``megacc`` for the MEGA-CC
    bridge) invoke whatever release the user has installed: no version is pinned
    anywhere in the code or CI, so the only way a result can be attributed to an
    inference program is to record what that program says about itself.  The
    probe is advisory: ``check=False``, non-zero exits and empty output are
    reported in the returned notes rather than raised, and the returned string
    says outright when nothing could be read, so a report can never be mistaken
    for "version recorded".

    ``runner`` defaults to :func:`run_tool`; callers pass their module-level
    binding so a test stub installed on that module intercepts the probe too.
    """
    run = runner or run_tool
    try:
        probe = run([executable, "--version"], check=False)
    except (OSError, ValueError) as exc:  # pragma: no cover - defensive
        return (
            f"not recorded ({tool} could not be probed)",
            [f"could not probe the {tool} version: {exc}"],
        )
    stdout = getattr(probe, "stdout", b"") or b""
    stderr = getattr(probe, "stderr", b"") or b""
    raw = (stdout or stderr).decode(errors="replace").strip()
    code = getattr(probe, "returncode", "?")
    if not raw:
        return (
            f"not recorded (empty {tool} --version output)",
            [
                f"the {tool} --version probe returned no output (exit code "
                f"{code}); the release used is not recorded"
            ],
        )
    first = raw.splitlines()[0].strip()
    notes: list[str] = []
    if code:
        notes.append(
            f"{tool} --version exited {code}; banner read as best effort: "
            f"{first!r}"
        )
    return first, notes
