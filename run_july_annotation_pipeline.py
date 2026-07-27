#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/ukb/run_july_annotation_pipeline.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "run_july_annotation_pipeline.py moved to scripts/ukb/run_july_annotation_pipeline.py",
        DeprecationWarning,
        stacklevel=1,
    )
    target = (
        Path(__file__).resolve().parent / "scripts" / "ukb" / "run_july_annotation_pipeline.py"
    )
    runpy.run_path(str(target), run_name="__main__")
