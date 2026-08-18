#!/usr/bin/env python3
"""Thin CLI for GSE87571 epigenetic-clock validation.

Prefer: ``uv run rogen-clock evaluate-gse87571``
"""

from __future__ import annotations

import typer

from rogen_aging.clock.gse87571 import main

if __name__ == "__main__":
    typer.run(main)
