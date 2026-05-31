#!/usr/bin/env python3
"""Thin CLI for epigenetic clock training and evaluation."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from rogen_aging.clock.evaluate import evaluate_clock
from rogen_aging.clock.train import train_clock


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Train or evaluate an epigenetic clock.")
    sub = parser.add_subparsers(dest="command", required=True)

    train_p = sub.add_parser("train", help="Fit ElasticNetCV clock on a wide methylation table.")
    train_p.add_argument("--input_data", type=Path, required=True, help="Parquet/CSV with cg* + chronological_age.")
    train_p.add_argument("--output_model", type=Path, required=True, help="Path for the fitted pipeline (.pkl/.joblib).")
    train_p.add_argument("--output_metrics", type=Path, required=True, help="Path for training metrics JSON.")
    train_p.add_argument("--test_size", type=float, default=0.2, help="Held-out test fraction (default 0.2).")
    train_p.add_argument("--random_state", type=int, default=42, help="Random seed (default 42).")

    eval_p = sub.add_parser("evaluate", help="Validate a saved clock on held-out data.")
    eval_p.add_argument("--model_path", type=Path, required=True, help="Trained model (.pkl or .joblib).")
    eval_p.add_argument("--test_data", type=Path, required=True, help="Test table (.parquet or .csv).")
    eval_p.add_argument("--output_dir", type=Path, required=True, help="Directory for figures and metrics JSON.")

    return parser


def main(argv: list[str] | None = None) -> None:
    """Parse CLI arguments and dispatch to :mod:`rogen_aging.clock`."""
    args = _build_parser().parse_args(argv)
    if args.command == "train":
        logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
        metrics = train_clock(
            args.input_data,
            args.output_model,
            args.output_metrics,
            test_size=args.test_size,
            random_state=args.random_state,
        )
        print(
            f"# CpGs used: {metrics['n_cpgs_features']} | alpha: {metrics['alpha']:.6g} | "
            f"l1_ratio: {metrics['l1_ratio']:.4g} | test MAE: {metrics['test_mae']:.4f} | "
            f"test r: {metrics['test_pearson_r']:.4f}"
        )
        return

    if args.command == "evaluate":
        result = evaluate_clock(args.model_path, args.test_data, args.output_dir)
        imputed = result.pop("imputed_missing_cpgs", [])
        print(json.dumps(result, indent=2))
        if imputed:
            print(f"Imputed {len(imputed)} missing model CpGs (see metrics JSON for names).", file=sys.stderr)
        return

    raise SystemExit(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
