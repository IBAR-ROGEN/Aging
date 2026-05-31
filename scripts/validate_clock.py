#!/usr/bin/env python3
"""Validate a pre-trained epigenetic clock (e.g. ElasticNet) on held-out bedMethyl-style data.

Backward-compatible CLI wrapper around :mod:`rogen_aging.clock.evaluate`.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from rogen_aging.clock.evaluate import evaluate_clock


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate a pre-trained epigenetic clock on held-out methylation + age data.",
    )
    parser.add_argument(
        "--model_path",
        type=Path,
        required=True,
        help="Path to a trained model (.pkl or .joblib).",
    )
    parser.add_argument(
        "--test_data",
        type=Path,
        required=True,
        help="Path to test table (.parquet or .csv) with chronological_age and cg* CpG columns.",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        required=True,
        help="Directory for figures and metrics JSON (created if missing).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = evaluate_clock(args.model_path, args.test_data, args.output_dir)
    imputed = result.pop("imputed_missing_cpgs", [])
    print(json.dumps({k: v for k, v in result.items() if k != "imputed_missing_cpgs"}, indent=2))
    if imputed:
        print(f"Imputed {len(imputed)} missing model CpGs (see metrics JSON for names).", file=sys.stderr)


if __name__ == "__main__":
    main()
