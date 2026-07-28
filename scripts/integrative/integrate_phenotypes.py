#!/usr/bin/env python3
"""Compute composite phenotypic risk from molecularly annotated variants.

Production default reads ``analysis/integrative/results/annotated_variants.parquet``
(from the tissue-map / full pipeline) and writes risk profiles beside it.

Example:
    uv run python scripts/integrative/integrate_phenotypes.py
    uv run python scripts/integrative/integrate_phenotypes.py --demo
"""

from __future__ import annotations

from pathlib import Path

import typer

from rogen_aging.integrative import PhenotypeIntegrator, VariantTissueMapper
from rogen_aging.integrative.io import (
    DEFAULT_OUTPUT_DIR,
    REPO_ROOT,
    read_table,
)

app = typer.Typer(add_completion=False, help=__doc__)

DEFAULT_ANNOTATED = DEFAULT_OUTPUT_DIR / "annotated_variants.parquet"


@app.command()
def main(
    annotated: Path | None = typer.Option(
        None,
        "--annotated",
        help=(
            "Tissue-mapped annotated variant table. "
            f"Default: {DEFAULT_ANNOTATED.relative_to(REPO_ROOT)}."
        ),
    ),
    output_dir: Path = typer.Option(
        DEFAULT_OUTPUT_DIR,
        "--output-dir",
        "-o",
        help="Directory for risk-profile Parquet outputs.",
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
    """Link molecular scores with composite phenotypic risk profiles.

    Args:
        annotated: Path to tissue-mapped annotated variants (Parquet/CSV).
        output_dir: Directory for ``variant_risks`` / ``sample_profiles`` Parquet.
        samples: Optional long genotype table for sample-level aggregation.
        demo: Materialize fixtures and score them.
    """
    if demo:
        from rogen_aging.pipeline_fixtures import write_integrative_fixtures

        fixtures = write_integrative_fixtures(repo_root=REPO_ROOT)
        mapped = VariantTissueMapper().map_variants_to_tissues(
            read_table(fixtures["variants"]),
            read_table(fixtures["eqtls"]),
            probe_annotation=read_table(fixtures["probes"]),
        )
        annotated_path = (
            REPO_ROOT / "analysis" / "integrative" / "fixtures" / "mapped_variants.parquet"
        )
        mapped["annotated"].write_parquet(annotated_path)
        annotated = annotated_path
        samples = fixtures["samples"]
        if output_dir.resolve() == DEFAULT_OUTPUT_DIR.resolve():
            output_dir = REPO_ROOT / "analysis" / "integrative" / "demo"
    elif annotated is None:
        annotated = DEFAULT_ANNOTATED
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
    output_dir.mkdir(parents=True, exist_ok=True)
    profiles["variant_risks"].write_parquet(output_dir / "variant_risks.parquet")
    if "sample_profiles" in profiles:
        profiles["sample_profiles"].write_parquet(output_dir / "sample_profiles.parquet")
    typer.echo(
        f"Wrote phenotype risk profiles | variants={profiles['variant_risks'].height} "
        f"| output={output_dir.resolve()}"
    )


if __name__ == "__main__":
    app()
