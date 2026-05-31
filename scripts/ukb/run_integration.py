#!/usr/bin/env python3
"""Activity 2.1.11.1 — thin CLI for synthetic UKB integrative validation.

Joins mock RAP phenotype CSV and LA-SNP VCF (from ``rogen-ukb-mock-rap``), runs
dominant-model association scans, and writes tidy result tables. **Synthetic data
only**; no biological conclusions.

Prefer: ``uv run rogen-ukb-integrate``
"""

from __future__ import annotations

from rogen_aging.integration.run_cli import main

if __name__ == "__main__":
    raise SystemExit(main())
