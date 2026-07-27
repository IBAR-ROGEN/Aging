#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/vcf/generate_synthetic_romanian_vcf.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "scripts/generate_synthetic_romanian_vcf.py moved; prefer `uv run rogen-vcf-synthetic`.",
        DeprecationWarning,
        stacklevel=1,
    )
    runpy.run_path(
        str(Path(__file__).resolve().parent / "vcf" / "generate_synthetic_romanian_vcf.py"),
        run_name="__main__",
    )
