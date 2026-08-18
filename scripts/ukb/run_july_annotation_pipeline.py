#!/usr/bin/env python3
"""Thin CLI for July prioritized-variant annotation.

Prefer: ``uv run rogen-july-annotate``
"""

from __future__ import annotations

from rogen_aging.annotation.july import app

if __name__ == "__main__":
    app()
