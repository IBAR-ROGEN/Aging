#!/usr/bin/env python3
"""Train an ElasticNetCV epigenetic clock from a sample-by-CpG table plus chronological age.

GSE40279 (Hannum 2013, Illumina 450K whole blood) must be obtained from GEO and converted
to a table with rows=samples, CpG columns whose names start with ``cg``, and a
``chronological_age`` column. That download and conversion are not implemented here.

Implementation lives in :mod:`rogen_aging.clock`; this script is a backward-compatible CLI wrapper.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from rogen_aging.clock.train import train_clock


def parse_args() -> argparse.Namespace:
    """CLI: input table, output paths, split and seed."""
    p = argparse.ArgumentParser(description="Train ElasticNetCV clock on methylation + age table.")
    p.add_argument("--input_data", type=Path, required=True, help="Parquet or CSV: samples x (cg* + chronological_age).")
    p.add_argument("--output_model", type=Path, required=True, help="Path to write trained model (.pkl).")
    p.add_argument("--output_metrics", type=Path, required=True, help="Path to write training metrics JSON.")
    p.add_argument("--test_size", type=float, default=0.2, help="Fraction held out for test (default 0.2).")
    p.add_argument("--random_state", type=int, default=42, help="Random seed (default 42).")
    return p.parse_args()


def main() -> None:
    """Train/test split, fit pipeline, write model + metrics JSON, print summary."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = parse_args()
    metrics = train_clock(
        args.input_data,
        args.output_model,
        args.output_metrics,
        test_size=args.test_size,
        random_state=args.random_state,
    )
    print(
        f"# CpGs used: {metrics['n_cpgs_features']} | alpha: {metrics['alpha']:.6g} | l1_ratio: {metrics['l1_ratio']:.4g} | "
        f"test MAE: {metrics['test_mae']:.4f} | test r: {metrics['test_pearson_r']:.4f}"
    )


if __name__ == "__main__":
    main()
