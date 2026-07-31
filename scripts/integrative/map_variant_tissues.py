#!/usr/bin/env python3
"""Offline variant→tissue mapping using ``rogen_aging.integrative``.

Production default reads the July annotation Combined_Master + GTEx long table
and writes Parquet outputs under ``analysis/integrative/results/``.

Example:
    uv run python scripts/integrative/map_variant_tissues.py
    uv run python scripts/integrative/map_variant_tissues.py --demo
    uv run python scripts/integrative/map_variant_tissues.py --config config/production.yaml
"""

from __future__ import annotations

from pathlib import Path

import typer

from rogen_aging.config import cfg_path, find_repo_root
from rogen_aging.config.cli import config_option, load_cli_config
from rogen_aging.integrative import VariantTissueMapper
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
        help="Directory for annotated + summary Parquet outputs. Default: from config.",
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
    """Map annotated variants onto tissue eQTL (+ optional methylation) profiles."""
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
        if output_dir is None or resolved_output.resolve() == default_output.resolve():
            resolved_output = demo_output
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
    resolved_output.mkdir(parents=True, exist_ok=True)
    result["annotated"].write_parquet(resolved_output / "annotated_variants.parquet")
    result["eqtl_summary"].write_parquet(resolved_output / "eqtl_summary.parquet")
    if "methylation_links" in result:
        result["methylation_links"].write_parquet(resolved_output / "methylation_links.parquet")
    typer.echo(
        f"Wrote integrative tissue map | variants={result['annotated'].height} "
        f"| output={resolved_output.resolve()}"
    )


if __name__ == "__main__":
    app()
