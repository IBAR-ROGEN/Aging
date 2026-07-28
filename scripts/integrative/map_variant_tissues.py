#!/usr/bin/env python3
"""Offline variant→tissue mapping using ``rogen_aging.integrative``.

Production default reads the July annotation Combined_Master + GTEx long table
and writes Parquet outputs under ``analysis/integrative/results/``.

Example:
    uv run python scripts/integrative/map_variant_tissues.py
    uv run python scripts/integrative/map_variant_tissues.py --demo
"""

from __future__ import annotations

from pathlib import Path

import typer

from rogen_aging.integrative import VariantTissueMapper
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
        help="Directory for annotated + summary Parquet outputs.",
    ),
    alphagenome: Path | None = typer.Option(
        None, "--alphagenome", help="Optional AlphaGenome score matrix."
    ),
    probes: Path | None = typer.Option(
        None, "--probes", help="Optional HM450/EPIC probe annotation CSV."
    ),
    demo: bool = typer.Option(
        False,
        "--demo",
        help="Write offline fixtures and map variants against them.",
    ),
) -> None:
    """Map annotated variants onto tissue eQTL (+ optional methylation) profiles.

    Args:
        variants: Path to the VEP/AlphaGenome annotated variant table.
        eqtls: Path to the long GTEx eQTL hit table.
        output_dir: Directory for ``annotated_variants`` / ``eqtl_summary`` Parquet.
        alphagenome: Optional AlphaGenome score matrix path.
        probes: Optional HM450/EPIC probe annotation path.
        demo: Materialize fixtures and run offline.
    """
    if demo:
        from rogen_aging.pipeline_fixtures import write_integrative_fixtures

        fixtures = write_integrative_fixtures(repo_root=REPO_ROOT)
        variants = fixtures["variants"]
        eqtls = fixtures["eqtls"]
        probes = fixtures["probes"]
        if output_dir.resolve() == DEFAULT_OUTPUT_DIR.resolve():
            output_dir = REPO_ROOT / "analysis" / "integrative" / "demo"
        variant_df = read_table(variants)
        eqtl_df = read_table(eqtls)
        alphagenome_df = None
        probe_df = read_table(probes) if probes else None
    else:
        ensure_july_parquet_cache()
        variant_df = load_production_variants(variants)
        eqtl_df = load_production_eqtls(eqtls)
        alphagenome_df = read_table(alphagenome) if alphagenome else None
        probe_df = read_table(probes) if probes else None

    mapper = VariantTissueMapper()
    result = mapper.map_variants_to_tissues(
        variant_df,
        eqtl_df,
        alphagenome=alphagenome_df,
        probe_annotation=probe_df,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    result["annotated"].write_parquet(output_dir / "annotated_variants.parquet")
    result["eqtl_summary"].write_parquet(output_dir / "eqtl_summary.parquet")
    if "methylation_links" in result:
        result["methylation_links"].write_parquet(output_dir / "methylation_links.parquet")
    typer.echo(
        f"Wrote integrative tissue map | variants={result['annotated'].height} "
        f"| output={output_dir.resolve()}"
    )


if __name__ == "__main__":
    app()
