#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/ukb/run_integration.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "scripts/run_integration.py moved to scripts/ukb/run_integration.py; "
        "prefer `uv run rogen-ukb-integrate`.",
        DeprecationWarning,
        stacklevel=1,
    )
    runpy.run_path(str(Path(__file__).resolve().parent / "ukb" / "run_integration.py"), run_name="__main__")
