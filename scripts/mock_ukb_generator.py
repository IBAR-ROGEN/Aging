#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/ukb/mock_clinical_csv.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "scripts/mock_ukb_generator.py moved to scripts/ukb/mock_clinical_csv.py; "
        "prefer `uv run rogen-ukb-mock-clinical`.",
        DeprecationWarning,
        stacklevel=1,
    )
    runpy.run_path(str(Path(__file__).resolve().parent / "ukb" / "mock_clinical_csv.py"), run_name="__main__")
