#!/usr/bin/env python3
"""Offline variant→tissue mapping using ``rogen_aging.integrative``.

Reads annotated variants and a long GTEx eQTL table, optionally joins
AlphaGenome scores and methylation probe annotations, and writes Parquet
outputs. Network calls stay in the legacy ``scripts/ukb/annotate_la_snps_*.py``
fetchers; this script only performs the integrative joins.

Example:
    uv run python scripts/integrative/map_variant_tissues.py \\
        --variants data/processed/prioritized_variants.csv \\
        --eqtls analysis/gtex_annotation/la_snp_gtex_eqtls.csv \\
        --output-dir analysis/integrative/
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import typer

from rogen_aging.integrative import VariantTissueMapper

app = typer.Typer(add_completion=False, help=__doc__)


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
    variants: Path = typer.Option(..., "--variants", help="Annotated variant table."),
    eqtls: Path = typer.Option(..., "--eqtls", help="Long GTEx eQTL CSV/Parquet."),
    output_dir: Path = typer.Option(
        Path("analysis/integrative"),
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
) -> None:
    """Map annotated variants onto tissue eQTL (+ optional methylation) profiles.

    Args:
        variants: Path to the VEP/AlphaGenome annotated variant table.
        eqtls: Path to the long GTEx eQTL hit table.
        output_dir: Directory for ``annotated_variants`` / ``eqtl_summary`` Parquet.
        alphagenome: Optional AlphaGenome score matrix path.
        probes: Optional HM450/EPIC probe annotation path.
    """
    mapper = VariantTissueMapper()
    result = mapper.map_variants_to_tissues(
        _read_table(variants),
        _read_table(eqtls),
        alphagenome=_read_table(alphagenome) if alphagenome else None,
        probe_annotation=_read_table(probes) if probes else None,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    result["annotated"].write_parquet(output_dir / "annotated_variants.parquet")
    result["eqtl_summary"].write_parquet(output_dir / "eqtl_summary.parquet")
    if "methylation_links" in result:
        result["methylation_links"].write_parquet(output_dir / "methylation_links.parquet")
    typer.echo(f"Wrote integrative tissue map under {output_dir.resolve()}")


if __name__ == "__main__":
    app()
