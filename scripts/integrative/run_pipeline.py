#!/usr/bin/env python3
"""End-to-end integrative pipeline CLI (variant→tissue→phenotype risk).

Production default (Activity A.2.1.11.1) reads the July annotation workbook /
parquet and writes results under ``analysis/integrative/results/``.

Example:
    uv run python scripts/integrative/run_pipeline.py
    uv run python scripts/integrative/run_pipeline.py --demo
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import typer

from rogen_aging.integrative import run_integrative_pipeline
from rogen_aging.integrative.io import (
    DEFAULT_OUTPUT_DIR,
    REPO_ROOT,
    ensure_july_parquet_cache,
    load_production_eqtls,
    load_production_variants,
    read_table,
)

app = typer.Typer(add_completion=False, help=__doc__)


@app.command()
def main(
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
    output_dir: Path = typer.Option(
        DEFAULT_OUTPUT_DIR,
        "--output-dir",
        "-o",
        help="Output directory for pipeline Parquet artefacts.",
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

    Args:
        variants: Path to the annotated variant table (production default if omitted).
        eqtls: Path to the long GTEx eQTL table (production default if omitted).
        output_dir: Directory for output Parquet files.
        alphagenome: Optional AlphaGenome score matrix.
        probes: Optional HM450/EPIC probe→gene annotation.
        samples: Optional long genotype table for sample-level risk.
        demo: Materialize fixtures and run offline.
    """
    if demo:
        from rogen_aging.pipeline_fixtures import write_integrative_fixtures

        fixtures = write_integrative_fixtures(repo_root=REPO_ROOT)
        variants = fixtures["variants"]
        eqtls = fixtures["eqtls"]
        probes = fixtures["probes"]
        samples = fixtures["samples"]
        if output_dir.resolve() == DEFAULT_OUTPUT_DIR.resolve():
            output_dir = REPO_ROOT / "analysis" / "integrative" / "demo"
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
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    for key in ("annotated", "eqtl_summary", "variant_risks", "methylation_links", "sample_profiles"):
        frame = result.get(key)
        if isinstance(frame, pl.DataFrame):
            out_name = "annotated_variants.parquet" if key == "annotated" else f"{key}.parquet"
            frame.write_parquet(output_dir / out_name)
            written.append(out_name)
    typer.echo(
        f"Integrative pipeline complete | variants={result['variant_risks'].height} "
        f"| files={','.join(written)} | output={output_dir.resolve()}"
    )


if __name__ == "__main__":
    app()
