"""``rogen-clock`` console entry — train and evaluate epigenetic clocks."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import typer

from rogen_aging.clock.evaluate import evaluate_clock
from rogen_aging.clock.train import train_clock
from rogen_aging.config import get_config
from rogen_aging.config.cli import config_option, load_cli_config

app = typer.Typer(
    add_completion=False, no_args_is_help=True, help="Train or evaluate an epigenetic clock."
)


@app.command("train")
def train_cmd(
    input_data: Path = typer.Option(
        ..., "--input_data", help="Parquet/CSV with cg* + chronological_age."
    ),
    output_model: Path = typer.Option(
        ..., "--output_model", help="Path for fitted pipeline (.pkl/.joblib)."
    ),
    output_metrics: Path = typer.Option(
        ..., "--output_metrics", help="Training metrics JSON path."
    ),
    config: Path | None = config_option(),
    test_size: float | None = typer.Option(
        None, "--test_size", help="Held-out test fraction. Default: from config."
    ),
    random_state: int | None = typer.Option(
        None, "--random_state", help="Random seed. Default: from config."
    ),
) -> None:
    """Train an epigenetic clock and write model plus metrics."""
    load_cli_config(config)
    cfg = get_config().clock
    resolved_test_size = float(cfg.test_size if test_size is None else test_size)
    resolved_random_state = int(cfg.random_state if random_state is None else random_state)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    metrics = train_clock(
        input_data,
        output_model,
        output_metrics,
        test_size=resolved_test_size,
        random_state=resolved_random_state,
    )
    typer.echo(
        f"# CpGs used: {metrics['n_cpgs_features']} | alpha: {metrics['alpha']:.6g} | "
        f"l1_ratio: {metrics['l1_ratio']:.4g} | test MAE: {metrics['test_mae']:.4f} | "
        f"test r: {metrics['test_pearson_r']:.4f}"
    )


@app.command("evaluate")
def evaluate_cmd(
    model_path: Path = typer.Option(..., "--model_path", help="Trained model (.pkl or .joblib)."),
    test_data: Path = typer.Option(..., "--test_data", help="Test table (.parquet or .csv)."),
    output_dir: Path = typer.Option(
        ..., "--output_dir", help="Directory for figures and metrics JSON."
    ),
    config: Path | None = config_option(),
) -> None:
    """Evaluate a trained clock on held-out data and write figures/metrics."""
    load_cli_config(config)
    result = evaluate_clock(model_path, test_data, output_dir)
    imputed = result.pop("imputed_missing_cpgs", [])
    typer.echo(json.dumps(result, indent=2))
    if imputed:
        typer.echo(
            f"Imputed {len(imputed)} missing model CpGs (see metrics JSON for names).",
            err=True,
        )


def main() -> None:
    """Console entry for ``rogen-clock``."""
    app()


if __name__ == "__main__":
    main()
