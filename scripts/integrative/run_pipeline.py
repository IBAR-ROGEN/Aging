#!/usr/bin/env python3
"""End-to-end integrative pipeline CLI (variant→tissue→phenotype risk).

Production default (Activity A.2.1.11.1) reads the July annotation workbook /
parquet and writes results under ``analysis/integrative/results/``.

Example:
    uv run python scripts/integrative/run_pipeline.py
    uv run python scripts/integrative/run_pipeline.py --demo
    uv run python scripts/integrative/run_pipeline.py --config config/production.yaml
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import typer

from rogen_aging.config import cfg_path, find_repo_root
from rogen_aging.config.cli import config_option, load_cli_config
from rogen_aging.integrative import run_integrative_pipeline
from rogen_aging.integrative.io import (
    ensure_july_parquet_cache,
    load_production_eqtls,
    load_production_variants,
    read_table,
)

app = typer.Typer(add_completion=False, help=__doc__)


@app.command()
def main(
    config: Path | None = config_option(),
    variants: Path | None = typer.Option(
        None,
        "--variants",
        help=(
            "Annotated variant table. Default: July Combined_Master parquet "
            "(or Supplementary_Table_1_Annotated_Variants.xlsx)."
        ),
    ),
    eqtls: Path | None = typer.Option(
        None,
        "--eqtls",
        help=(
            "Long GTEx eQTL table. Default: July GTEx parquet "
            "(or analysis/gtex_annotation/la_snp_gtex_eqtls.csv)."
        ),
    ),
    output_dir: Path | None = typer.Option(
        None,
        "--output-dir",
        "-o",
        help="Output directory for pipeline Parquet artefacts. Default: from config.",
    ),
    alphagenome: Path | None = typer.Option(
        None,
        "--alphagenome",
        help="Optional AlphaGenome scores (usually already on Combined_Master).",
    ),
    probes: Path | None = typer.Option(None, "--probes"),
    samples: Path | None = typer.Option(None, "--samples"),
    demo: bool = typer.Option(
        False,
        "--demo",
        help="Write offline fixtures and run the pipeline against them.",
    ),
) -> None:
    """Execute the full integrative variant→tissue→phenotype pipeline.

    Without ``--demo``, defaults to the July production annotation artefacts and
    writes Parquet outputs under ``analysis/integrative/results/``.
    """
    cfg = load_cli_config(config)
    repo_root = find_repo_root()
    default_output = cfg_path(cfg, "paths", "integrative", "output_dir")
    demo_output = cfg_path(cfg, "paths", "integrative", "demo_dir")
    resolved_output = output_dir or default_output

    if demo:
        from rogen_aging.pipeline_fixtures import write_integrative_fixtures

        fixtures = write_integrative_fixtures(repo_root=repo_root)
        variants = fixtures["variants"]
        eqtls = fixtures["eqtls"]
        probes = fixtures["probes"]
        samples = fixtures["samples"]
        if output_dir is None or resolved_output.resolve() == default_output.resolve():
            resolved_output = demo_output
        variant_df = read_table(variants)
        eqtl_df = read_table(eqtls)
        alphagenome_df = None
        probe_df = read_table(probes) if probes else None
        sample_df = read_table(samples) if samples else None
    else:
        ensure_july_parquet_cache()
        variant_df = load_production_variants(variants)
        eqtl_df = load_production_eqtls(eqtls)
        alphagenome_df = read_table(alphagenome) if alphagenome else None
        probe_df = read_table(probes) if probes else None
        sample_df = read_table(samples) if samples else None

    result = run_integrative_pipeline(
        variant_df,
        eqtl_df,
        alphagenome=alphagenome_df,
        probe_annotation=probe_df,
        sample_phenotypes=sample_df,
    )
    resolved_output.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    for key in (
        "annotated",
        "eqtl_summary",
        "variant_risks",
        "methylation_links",
        "sample_profiles",
    ):
        frame = result.get(key)
        if isinstance(frame, pl.DataFrame):
            out_name = "annotated_variants.parquet" if key == "annotated" else f"{key}.parquet"
            frame.write_parquet(resolved_output / out_name)
            written.append(out_name)
    typer.echo(
        f"Integrative pipeline complete | variants={result['variant_risks'].height} "
        f"| files={','.join(written)} | output={resolved_output.resolve()}"
    )


if __name__ == "__main__":
    app()
