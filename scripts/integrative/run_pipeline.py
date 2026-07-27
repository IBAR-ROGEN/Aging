#!/usr/bin/env python3
"""End-to-end integrative pipeline CLI (variant→tissue→phenotype risk).

Example:
    uv run python scripts/integrative/run_pipeline.py --demo
    uv run python scripts/integrative/run_pipeline.py \\
        --variants analysis/integrative/fixtures/annotated_variants.parquet \\
        --eqtls analysis/integrative/fixtures/eqtls.parquet \\
        --output-dir analysis/integrative/
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import typer

from rogen_aging.integrative import run_integrative_pipeline

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
    variants: Path | None = typer.Option(
        None,
        "--variants",
        help="Annotated variant table (required unless --demo).",
    ),
    eqtls: Path | None = typer.Option(
        None,
        "--eqtls",
        help="Long GTEx eQTL table (required unless --demo).",
    ),
    output_dir: Path = typer.Option(
        Path("analysis/integrative"),
        "--output-dir",
        "-o",
        help="Output directory for pipeline Parquet artefacts.",
    ),
    alphagenome: Path | None = typer.Option(None, "--alphagenome"),
    probes: Path | None = typer.Option(None, "--probes"),
    samples: Path | None = typer.Option(None, "--samples"),
    demo: bool = typer.Option(
        False,
        "--demo",
        help="Write offline fixtures and run the pipeline against them.",
    ),
) -> None:
    """Execute the full integrative variant→tissue→phenotype pipeline.

    Writes Parquet artefacts (``annotated``, ``eqtl_summary``, ``variant_risks``,
    and optionally ``methylation_links`` / ``sample_profiles``) under
    ``output_dir``.

    Args:
        variants: Path to the annotated variant table.
        eqtls: Path to the long GTEx eQTL table.
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
        if output_dir == Path("analysis/integrative"):
            output_dir = REPO_ROOT / "analysis" / "integrative" / "demo"

    if variants is None or eqtls is None:
        raise typer.BadParameter("Provide --variants and --eqtls, or pass --demo.")

    result = run_integrative_pipeline(
        _read_table(variants),
        _read_table(eqtls),
        alphagenome=_read_table(alphagenome) if alphagenome else None,
        probe_annotation=_read_table(probes) if probes else None,
        sample_phenotypes=_read_table(samples) if samples else None,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    for key in ("annotated", "eqtl_summary", "variant_risks", "methylation_links", "sample_profiles"):
        frame = result.get(key)
        if isinstance(frame, pl.DataFrame):
            frame.write_parquet(output_dir / f"{key}.parquet")
    typer.echo(
        f"Integrative pipeline complete | variants={result['variant_risks'].height} "
        f"| output={output_dir.resolve()}"
    )


if __name__ == "__main__":
    app()
