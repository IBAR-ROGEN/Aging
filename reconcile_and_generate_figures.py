#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/figures/reconcile_and_generate_figures.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "reconcile_and_generate_figures.py moved to "
        "scripts/figures/reconcile_and_generate_figures.py",
        DeprecationWarning,
        stacklevel=1,
    )
    target = (
        Path(__file__).resolve().parent
        / "scripts"
        / "figures"
        / "reconcile_and_generate_figures.py"
    )
    runpy.run_path(str(target), run_name="__main__")
