#!/usr/bin/env python3
"""Compute composite phenotypic risk from molecularly annotated variants.

Example:
    uv run python scripts/integrative/integrate_phenotypes.py \\
        --annotated analysis/integrative/annotated_variants.parquet \\
        --output-dir analysis/integrative/
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import typer

from rogen_aging.integrative import PhenotypeIntegrator

app = typer.Typer(add_completion=False, help=__doc__)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent


def _read_table(path: Path) -> pl.DataFrame:
    """Load a CSV, TSV, Excel, or Parquet table into Polars.

    Args:
        path: Input table path. Format is inferred from the file suffix
            (``.parquet``, ``.xlsx``/``.xls``, ``.tsv``, otherwise CSV).

    Returns:
        Polars DataFrame containing the loaded table.
    """
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pl.read_parquet(path)
    if suffix in {".xlsx", ".xls"}:
        return pl.from_pandas(__import__("pandas").read_excel(path))
    if suffix == ".tsv":
        return pl.read_csv(path, separator="\t")
    return pl.read_csv(path)


@app.command()
def main(
    annotated: Path | None = typer.Option(
        None, "--annotated", help="Tissue-mapped annotated variant table (or --demo)."
    ),
    output_dir: Path = typer.Option(
        Path("analysis/integrative"),
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
        from rogen_aging.integrative import VariantTissueMapper
        from rogen_aging.pipeline_fixtures import write_integrative_fixtures

        fixtures = write_integrative_fixtures(repo_root=REPO_ROOT)
        mapped = VariantTissueMapper().map_variants_to_tissues(
            _read_table(fixtures["variants"]),
            _read_table(fixtures["eqtls"]),
            probe_annotation=_read_table(fixtures["probes"]),
        )
        annotated_path = (
            REPO_ROOT / "analysis" / "integrative" / "fixtures" / "mapped_variants.parquet"
        )
        mapped["annotated"].write_parquet(annotated_path)
        annotated = annotated_path
        samples = fixtures["samples"]
        if output_dir == Path("analysis/integrative"):
            output_dir = REPO_ROOT / "analysis" / "integrative" / "demo"

    if annotated is None:
        raise typer.BadParameter("Provide --annotated, or pass --demo.")

    integrator = PhenotypeIntegrator()
    profiles = integrator.build_risk_profile(
        _read_table(annotated),
        sample_phenotypes=_read_table(samples) if samples else None,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    profiles["variant_risks"].write_parquet(output_dir / "variant_risks.parquet")
    if "sample_profiles" in profiles:
        profiles["sample_profiles"].write_parquet(output_dir / "sample_profiles.parquet")
    typer.echo(f"Wrote phenotype risk profiles under {output_dir.resolve()}")


if __name__ == "__main__":
    app()
