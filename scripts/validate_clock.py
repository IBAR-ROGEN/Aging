#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/clock/validate_clock.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "scripts/validate_clock.py moved; prefer `uv run rogen-clock evaluate`.",
        DeprecationWarning,
        stacklevel=1,
    )
    runpy.run_path(str(Path(__file__).resolve().parent / "clock" / "validate_clock.py"), run_name="__main__")
