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
:mod:`openreltime.megacc` (MEGA-CC).
"""

from __future__ import annotations

import re
import subprocess  # noqa: S404 - argv-list execution only, shell=False
from pathlib import Path
from shutil import which

#: Characters that never belong in a command token.  With ``shell=False``
#: they could not be interpolated anyway; the blacklist is defense in depth
#: against accidental misuse of :func:`run_tool` with assembled strings.
#: Whitespace and non-ASCII characters are deliberately allowed (paths).
_FORBIDDEN = re.compile(r"[;&|<>$`\\\"'\r\n\t]")


def validate_tokens(*tokens: str) -> None:
    """Raise :class:`ValueError` if any token contains a shell metacharacter."""
    for token in tokens:
        if not token:
            raise ValueError("refusing to run: empty command token")
        if _FORBIDDEN.search(token):
            raise ValueError(
                f"refusing to run: argument {token!r} contains shell "
                "metacharacters (one of ; & | < > $ ` \\ \" ' or a newline)"
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
