#!/usr/bin/env python3
"""Deprecated wrapper — use ``uv run rogen-clock evaluate`` or ``scripts/clock/run_clock.py evaluate``."""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import typer

from rogen_aging.clock.evaluate import evaluate_clock

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Evaluate a pre-trained epigenetic clock on held-out methylation + age data.",
)


@app.callback(invoke_without_command=True)
def main(
    model_path: Path = typer.Option(..., "--model_path", help="Path to trained model (.pkl or .joblib)."),
    test_data: Path = typer.Option(..., "--test_data", help="Test table with chronological_age and cg* columns."),
    output_dir: Path = typer.Option(..., "--output_dir", help="Directory for figures and metrics JSON."),
) -> None:
    """Backward-compatible held-out validation wrapper."""
    warnings.warn(
        "scripts/validate_clock.py is deprecated; use `uv run rogen-clock evaluate` "
        "or `scripts/clock/run_clock.py evaluate`.",
        DeprecationWarning,
        stacklevel=1,
    )
    result = evaluate_clock(model_path, test_data, output_dir)
    imputed = result.pop("imputed_missing_cpgs", [])
    typer.echo(json.dumps(result, indent=2))
    if imputed:
        typer.echo(
            f"Imputed {len(imputed)} missing model CpGs (see metrics JSON for names).",
            err=True,
        )


if __name__ == "__main__":
    app(prog_name="validate_clock.py")
