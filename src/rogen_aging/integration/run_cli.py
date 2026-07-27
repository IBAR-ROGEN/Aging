"""Deprecated alias — use :mod:`rogen_aging.ukb_integration.run_cli`."""

from __future__ import annotations

import warnings

from rogen_aging.ukb_integration.run_cli import app, main

warnings.warn(
    "rogen_aging.integration.run_cli is deprecated; use rogen_aging.ukb_integration.run_cli",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["app", "main"]
