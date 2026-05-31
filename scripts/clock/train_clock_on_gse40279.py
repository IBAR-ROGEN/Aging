#!/usr/bin/env python3
"""Deprecated wrapper — use ``uv run rogen-clock train`` or ``scripts/clock/run_clock.py train``."""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import typer

from rogen_aging.clock.train import train_clock

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.callback(invoke_without_command=True)
def main(
    input_data: Path = typer.Option(..., "--input_data", help="Parquet or CSV with cg* + chronological_age."),
    output_model: Path = typer.Option(..., "--output_model", help="Path for trained model (.pkl)."),
    output_metrics: Path = typer.Option(..., "--output_metrics", help="Training metrics JSON path."),
    test_size: float = typer.Option(0.2, "--test_size"),
    random_state: int = typer.Option(42, "--random_state"),
) -> None:
    warnings.warn(
        "scripts/train_clock_on_gse40279.py is deprecated; use `uv run rogen-clock train`.",
        DeprecationWarning,
        stacklevel=1,
    )
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    metrics = train_clock(
        input_data,
        output_model,
        output_metrics,
        test_size=test_size,
        random_state=random_state,
    )
    typer.echo(
        f"# CpGs used: {metrics['n_cpgs_features']} | alpha: {metrics['alpha']:.6g} | "
        f"l1_ratio: {metrics['l1_ratio']:.4g} | test MAE: {metrics['test_mae']:.4f} | "
        f"test r: {metrics['test_pearson_r']:.4f}"
    )


if __name__ == "__main__":
    app(prog_name="train_clock_on_gse40279.py")
