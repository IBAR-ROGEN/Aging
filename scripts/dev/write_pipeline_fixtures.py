#!/usr/bin/env python3
"""Write offline fixtures so ROGEN pipelines can run without live data/APIs.

Example:
    uv run python scripts/dev/write_pipeline_fixtures.py
    uv run python scripts/ukb/run_july_annotation_pipeline.py --demo
"""

from __future__ import annotations

from pathlib import Path

import typer

from rogen_aging.pipeline_fixtures import write_all_pipeline_fixtures

app = typer.Typer(add_completion=False, help=__doc__)


@app.command()
def main(
    repo_root: Path = typer.Option(
        Path("."),
        "--repo-root",
        help="Repository root (defaults to cwd).",
    ),
    july_limit: int = typer.Option(
        12,
        "--july-limit",
        help="Max variants for the July demo table.",
    ),
    force_synthetic_clock: bool = typer.Option(
        False,
        "--force-synthetic-clock",
        help="Overwrite methylation inputs with a tiny synthetic cohort.",
    ),
) -> None:
    """Materialize July, integrative, and clock development fixtures.

    Args:
        repo_root: Repository root directory.
        july_limit: Maximum prioritized variants for the July demo CSV.
        force_synthetic_clock: Replace GSE methylation inputs with synthetic data.
    """
    paths = write_all_pipeline_fixtures(
        repo_root=repo_root.resolve(),
        july_limit=july_limit,
        force_synthetic_clock=force_synthetic_clock,
    )
    typer.echo("Wrote pipeline fixtures:")
    for name, path in (
        ("july_variants", paths.july_variants),
        ("july_alphagenome", paths.july_alphagenome),
        ("july_alphamissense", paths.july_alphamissense),
        ("july_local_vep", paths.july_local_vep),
        ("july_cache_dir", paths.july_cache_dir),
        ("integrative_variants", paths.integrative_variants),
        ("integrative_eqtls", paths.integrative_eqtls),
        ("integrative_probes", paths.integrative_probes),
        ("integrative_samples", paths.integrative_samples),
        ("clock_model", paths.clock_model),
        ("clock_methylation", paths.clock_methylation),
        ("clock_meta", paths.clock_meta),
    ):
        typer.echo(f"  {name}: {path}")


if __name__ == "__main__":
    app()
