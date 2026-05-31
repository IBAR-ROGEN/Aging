#!/usr/bin/env python3
"""CLI for LA-SNP manifest build and 1KG allele-frequency extraction."""

from __future__ import annotations

import sys

from rogen_aging.ukb.manifest import main

if __name__ == "__main__":
    raise SystemExit(main())
