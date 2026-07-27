"""Pytest configuration for headless-safe matplotlib usage."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
