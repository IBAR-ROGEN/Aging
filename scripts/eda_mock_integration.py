#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/eda/eda_mock_integration.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "scripts/eda_mock_integration.py moved to scripts/eda/eda_mock_integration.py.",
        DeprecationWarning,
        stacklevel=1,
    )
    runpy.run_path(str(Path(__file__).resolve().parent / "eda" / "eda_mock_integration.py"), run_name="__main__")
