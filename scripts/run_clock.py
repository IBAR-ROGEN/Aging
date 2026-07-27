#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/clock/run_clock.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "scripts/run_clock.py moved to scripts/clock/run_clock.py; prefer `uv run rogen-clock`.",
        DeprecationWarning,
        stacklevel=1,
    )
    runpy.run_path(str(Path(__file__).resolve().parent / "clock" / "run_clock.py"), run_name="__main__")
