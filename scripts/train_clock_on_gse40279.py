#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/clock/train_clock_on_gse40279.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "scripts/train_clock_on_gse40279.py moved; prefer `uv run rogen-clock train`.",
        DeprecationWarning,
        stacklevel=1,
    )
    runpy.run_path(
        str(Path(__file__).resolve().parent / "clock" / "train_clock_on_gse40279.py"),
        run_name="__main__",
    )
