"""``rogen-clock`` console entry — train and evaluate epigenetic clocks."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import typer

from rogen_aging.clock.evaluate import evaluate_clock
from rogen_aging.clock.train import train_clock
from rogen_aging.clock.validate_matrix import (
    MissingValuePolicy,
    validate_methylation_matrix,
    write_validation_outputs,
)

app = typer.Typer(add_completion=False, no_args_is_help=True, help="Train or evaluate an epigenetic clock.")

@app.command("train")
def train_cmd(
    input_data: Path = typer.Option(..., "--input_data", help="Parquet/CSV with cg* + chronological_age."),
    output_model: Path = typer.Option(..., "--output_model", help="Path for fitted pipeline (.pkl/.joblib)."),
    output_metrics: Path = typer.Option(..., "--output_metrics", help="Training metrics JSON path."),
    test_size: float = typer.Option(0.2, "--test_size", help="Held-out test fraction."),
    random_state: int = typer.Option(42, "--random_state", help="Random seed."),
) -> None:
    """Train an epigenetic clock and write model plus metrics.

    Args:
        input_data: Parquet/CSV with ``cg*`` columns and ``chronological_age``.
        output_model: Destination path for the fitted pipeline (``.pkl``/``.joblib``).
        output_metrics: Destination path for training metrics JSON.
        test_size: Held-out test fraction.
        random_state: Random seed for the train/test split.
    """
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


@app.command("evaluate")
def evaluate_cmd(
    model_path: Path = typer.Option(..., "--model_path", help="Trained model (.pkl or .joblib)."),
    test_data: Path = typer.Option(..., "--test_data", help="Test table (.parquet or .csv)."),
    output_dir: Path = typer.Option(..., "--output_dir", help="Directory for figures and metrics JSON."),
) -> None:
    """Evaluate a trained clock on held-out data and write figures/metrics.

    Args:
        model_path: Trained model path (``.pkl`` or ``.joblib``).
        test_data: Test table path (``.parquet`` or ``.csv``).
        output_dir: Directory for figures and metrics JSON.
    """
    result = evaluate_clock(model_path, test_data, output_dir)
    imputed = result.pop("imputed_missing_cpgs", [])
    typer.echo(json.dumps(result, indent=2))
    if imputed:
        typer.echo(
            f"Imputed {len(imputed)} missing model CpGs (see metrics JSON for names).",
            err=True,
        )


@app.command("validate-matrix")
def validate_matrix_cmd(
    matrix: Path = typer.Option(..., "--matrix", help="Wide beta matrix (.parquet/.csv/.tsv)."),
    metadata: Path = typer.Option(..., "--metadata", help="Sample metadata manifest."),
    expected_cpgs: Path | None = typer.Option(
        None,
        "--expected-cpgs",
        help="Optional text/CSV list of expected CpG IDs (one per line).",
    ),
    missing_policy: MissingValuePolicy = typer.Option(
        MissingValuePolicy.REPORT,
        "--missing-policy",
        help="How to handle missing betas: report|fail|impute_column_mean|drop_sites.",
    ),
    sample_id_col: str | None = typer.Option(
        None,
        "--sample-id-col",
        help="Sample ID column in the matrix (auto-detected if omitted).",
    ),
    metadata_sample_id_col: str | None = typer.Option(
        None,
        "--metadata-sample-id-col",
        help="Sample ID column in metadata (auto-detected if omitted).",
    ),
    allow_inclusive_beta: bool = typer.Option(
        False,
        "--allow-inclusive-beta",
        help="Allow beta values of exactly 0 or 1 (default requires strict (0, 1)).",
    ),
    log_out: Path | None = typer.Option(
        None,
        "--log-out",
        help="Write the human-readable diagnostic log to this path.",
    ),
    report_json: Path | None = typer.Option(
        None,
        "--report-json",
        help="Write the machine-readable validation report JSON.",
    ),
    cleaned_matrix: Path | None = typer.Option(
        None,
        "--cleaned-matrix",
        help="Optional path for the cleaned matrix after missing-value handling.",
    ),
) -> None:
    """Pre-flight check a methylation matrix before clock inference/training.

    Args:
        matrix: Wide beta matrix path.
        metadata: Metadata manifest path.
        expected_cpgs: Optional expected CpG list.
        missing_policy: Missing-value handling policy.
        sample_id_col: Optional matrix sample ID column name.
        metadata_sample_id_col: Optional metadata sample ID column name.
        allow_inclusive_beta: If true, accept betas in ``[0, 1]`` instead of ``(0, 1)``.
        log_out: Optional diagnostic log path.
        report_json: Optional JSON report path.
        cleaned_matrix: Optional cleaned matrix output path.
    """
    report, cleaned = validate_methylation_matrix(
        matrix,
        metadata,
        expected_cpgs=expected_cpgs,
        sample_id_col=sample_id_col,
        metadata_sample_id_col=metadata_sample_id_col,
        missing_policy=missing_policy,
        require_strict_beta=not allow_inclusive_beta,
    )
    write_validation_outputs(
        report,
        cleaned,
        log_path=log_out,
        report_json_path=report_json,
        cleaned_matrix_path=cleaned_matrix,
    )
    typer.echo(report.format_log())
    if not report.passed:
        raise typer.Exit(code=1)


def main() -> None:
    """Console entry for ``rogen-clock``."""
    app()

if __name__ == "__main__":
    main()
