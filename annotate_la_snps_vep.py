#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/ukb/annotate_la_snps_vep.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "annotate_la_snps_vep.py moved to scripts/ukb/annotate_la_snps_vep.py",
        DeprecationWarning,
        stacklevel=1,
    )
    target = Path(__file__).resolve().parent / "scripts" / "ukb" / "annotate_la_snps_vep.py"
    runpy.run_path(str(target), run_name="__main__")
