#!/usr/bin/env python3
"""Summarize 1KG vs gnomAD v4 NFE allele-frequency comparison for reporting."""

from __future__ import annotations

import sys

from rogen_aging.ukb.gnomad import summarize_main

if __name__ == "__main__":
    raise SystemExit(summarize_main())
