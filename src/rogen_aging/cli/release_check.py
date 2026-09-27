"""``rogen-release-check`` console entry."""

from __future__ import annotations

from pathlib import Path

import typer
from loguru import logger

from rogen_aging.config.cli import config_option, load_cli_config
from rogen_aging.privacy.release_guard import (
    ReleaseGuardError,
    load_table,
    min_cell_count_from_config,
    run_release_guard,
    write_report,
    write_table,
)

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Check an aggregate table before it leaves the processing zone.",
)


@app.command("run")
def run_cmd(
    table: Path = typer.Option(
        ...,
        "--table",
        exists=True,
        dir_okay=False,
        help="Aggregate table (CSV, TSV, or Parquet) to check.",
    ),
    table_kind: str = typer.Option(
        ...,
        "--table-kind",
        help="allele_counts (AC/AN), strata (phenotype-group counts), or summary (identifier checks only).",
    ),
    report: Path = typer.Option(
        ...,
        "--report",
        help="JSON path for the PASS/FAIL report.",
    ),
    n_samples: int | None = typer.Option(
        None,
        "--n-samples",
        help="Fail if the table has this many rows (one row per sample).",
    ),
    write_sanitized: Path | None = typer.Option(
        None,
        "--write-sanitized",
        help="Write a suppressed copy only when checks 1 and 2 pass.",
    ),
    config: Path | None = config_option(),
) -> None:
    """Run the release guard and write the JSON report."""
    cfg = load_cli_config(config)
    try:
        min_cell_count = min_cell_count_from_config(cfg)
        frame = load_table(table)
        result = run_release_guard(
            frame,
            table_kind=table_kind,
            min_cell_count=min_cell_count,
            n_samples=n_samples,
        )
    except ReleaseGuardError as exc:
        logger.error("{message}", message=str(exc))
        raise typer.Exit(code=1) from exc

    sanitized_written = False
    if write_sanitized is not None and result.sanitized is not None and not result.release_blocked:
        write_table(result.sanitized, write_sanitized)
        sanitized_written = True
    payload = result.report(table=str(table), sanitized_written=sanitized_written)
    write_report(payload, report)
    for item in result.checks:
        typer.echo(f"{item.status} {item.check}: {item.reason}")
    typer.echo(f"Wrote release report: {report}")
    if result.release_blocked:
        raise typer.Exit(code=1)


def main() -> None:
    """Console entry for ``rogen-release-check``."""
    app()


if __name__ == "__main__":
    main()
