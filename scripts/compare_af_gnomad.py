#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/ukb/compare_af_gnomad.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "scripts/compare_af_gnomad.py moved to scripts/ukb/compare_af_gnomad.py; "
        "prefer `uv run rogen-compare-af-gnomad`.",
        DeprecationWarning,
        stacklevel=1,
    )
    runpy.run_path(str(Path(__file__).resolve().parent / "ukb" / "compare_af_gnomad.py"), run_name="__main__")
