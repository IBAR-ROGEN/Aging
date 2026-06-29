#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/figures/plot_af_comparison.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "plot_af_comparison.py moved to scripts/figures/plot_af_comparison.py",
        DeprecationWarning,
        stacklevel=1,
    )
    target = Path(__file__).resolve().parent / "scripts" / "figures" / "plot_af_comparison.py"
    runpy.run_path(str(target), run_name="__main__")
