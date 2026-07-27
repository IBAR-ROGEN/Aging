#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/ukb/la_snp_lookup.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "scripts/ukb_la_snp_lookup.py moved to scripts/ukb/la_snp_lookup.py; "
        "prefer `uv run rogen-ukb-manifest`.",
        DeprecationWarning,
        stacklevel=1,
    )
    runpy.run_path(str(Path(__file__).resolve().parent / "ukb" / "la_snp_lookup.py"), run_name="__main__")
