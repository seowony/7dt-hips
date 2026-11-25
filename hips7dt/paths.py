"""Working-directory convention: every data product lives under WORK.

WORK = $HIPS7DT_WORK, default /data8/skyteam/7DT (maintainer, 2026-11-25).
A relative path given on the command line (--out, --index, --in-root ...)
is taken relative to WORK, not to the current directory; an absolute path
is used as is. The code itself stays in the git checkout.
"""
import os
from pathlib import Path

WORK = Path(os.environ.get("HIPS7DT_WORK", "/data8/skyteam/7DT"))


def resolve(p) -> str:
    """Absolute paths unchanged; relative ones under WORK."""
    if p in (None, ""):
        return p
    p = Path(p).expanduser()
    return str(p if p.is_absolute() else WORK / p)


def selftest() -> bool:
    return (resolve("/abs/x") == "/abs/x"
            and resolve("index/a.csv") == str(WORK / "index/a.csv")
            and resolve("") == "")
