#!/usr/bin/env python3
"""Deprecated path shim — forwards to ``scripts/ukb/mock_rap_folder.py``."""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn(
        "scripts/ukb_mock_gen.py moved to scripts/ukb/mock_rap_folder.py; "
        "prefer `uv run rogen-ukb-mock-rap`.",
        DeprecationWarning,
        stacklevel=1,
    )
    runpy.run_path(str(Path(__file__).resolve().parent / "ukb" / "mock_rap_folder.py"), run_name="__main__")
