#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/clock/evaluate_methylation_clock.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "evaluate_methylation_clock.py moved to scripts/clock/evaluate_methylation_clock.py",
        DeprecationWarning,
        stacklevel=1,
    )
    target = Path(__file__).resolve().parent / "scripts" / "clock" / "evaluate_methylation_clock.py"
    runpy.run_path(str(target), run_name="__main__")
