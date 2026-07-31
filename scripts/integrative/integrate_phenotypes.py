#!/usr/bin/env python3
"""Compute composite phenotypic risk from molecularly annotated variants.

Production default reads ``analysis/integrative/results/annotated_variants.parquet``
(from the tissue-map / full pipeline) and writes risk profiles beside it.

Example:
    uv run python scripts/integrative/integrate_phenotypes.py
    uv run python scripts/integrative/integrate_phenotypes.py --demo
    uv run python scripts/integrative/integrate_phenotypes.py --config config/production.yaml
"""

from __future__ import annotations

from pathlib import Path

import typer

from rogen_aging.config import cfg_path, find_repo_root
from rogen_aging.config.cli import config_option, load_cli_config
from rogen_aging.integrative import PhenotypeIntegrator, VariantTissueMapper
from rogen_aging.integrative.io import read_table

app = typer.Typer(add_completion=False, help=__doc__)


@app.command()
def main(
    config: Path | None = config_option(),
    annotated: Path | None = typer.Option(
        None,
        "--annotated",
        help=(
            "Tissue-mapped annotated variant table. "
            "Default: <integrative.output_dir>/annotated_variants.parquet."
        ),
    ),
    output_dir: Path | None = typer.Option(
        None,
        "--output-dir",
        "-o",
        help="Directory for risk-profile Parquet outputs. Default: from config.",
    ),
    samples: Path | None = typer.Option(
        None,
        "--samples",
        help="Optional long genotype table (sample_id, rsid, alt_dosage).",
    ),
    demo: bool = typer.Option(
        False,
        "--demo",
        help="Run against integrative fixtures (writes fixtures if needed).",
    ),
) -> None:
    """Link molecular scores with composite phenotypic risk profiles."""
    cfg = load_cli_config(config)
    repo_root = find_repo_root()
    default_output = cfg_path(cfg, "paths", "integrative", "output_dir")
    demo_output = cfg_path(cfg, "paths", "integrative", "demo_dir")
    fixtures_dir = cfg_path(cfg, "paths", "integrative", "fixtures_dir")
    default_annotated = default_output / "annotated_variants.parquet"
    resolved_output = output_dir or default_output

    if demo:
        from rogen_aging.pipeline_fixtures import write_integrative_fixtures

        fixtures = write_integrative_fixtures(repo_root=repo_root)
        mapped = VariantTissueMapper().map_variants_to_tissues(
            read_table(fixtures["variants"]),
            read_table(fixtures["eqtls"]),
            probe_annotation=read_table(fixtures["probes"]),
        )
        annotated_path = fixtures_dir / "mapped_variants.parquet"
        mapped["annotated"].write_parquet(annotated_path)
        annotated = annotated_path
        samples = fixtures["samples"]
        if output_dir is None or resolved_output.resolve() == default_output.resolve():
            resolved_output = demo_output
    elif annotated is None:
        annotated = default_annotated
        if not annotated.is_file():
            raise typer.BadParameter(
                f"Production annotated table not found: {annotated}. "
                "Run scripts/integrative/run_pipeline.py or map_variant_tissues.py first, "
                "or pass --annotated / --demo."
            )

    integrator = PhenotypeIntegrator()
    profiles = integrator.build_risk_profile(
        read_table(annotated),
        sample_phenotypes=read_table(samples) if samples else None,
    )
    resolved_output.mkdir(parents=True, exist_ok=True)
    profiles["variant_risks"].write_parquet(resolved_output / "variant_risks.parquet")
    if "sample_profiles" in profiles:
        profiles["sample_profiles"].write_parquet(resolved_output / "sample_profiles.parquet")
    typer.echo(
        f"Wrote phenotype risk profiles | variants={profiles['variant_risks'].height} "
        f"| output={resolved_output.resolve()}"
    )


if __name__ == "__main__":
    app()
