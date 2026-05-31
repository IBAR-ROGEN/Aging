#!/usr/bin/env python3
"""CLI for synthetic UK Biobank-style clinical CSV generation."""

from __future__ import annotations

from rogen_aging.ukb.mock_clinical import app

if __name__ == "__main__":
    app()
