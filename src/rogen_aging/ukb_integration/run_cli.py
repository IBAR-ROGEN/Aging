"""CLI for Activity 2.1.11.1 synthetic UKB integrative validation."""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any

import typer

from rogen_aging.ukb_integration.ukb_joiner import (
    ACTIVITY_ID,
    SYNTHETIC_DISCLAIMER,
    run_integration_pipeline,
)

DEFAULT_PHENO = Path("test_data/mock_ukb_rap/phenotypes/ukb_phenotypes.csv")
DEFAULT_VCF = Path("test_data/mock_ukb_rap/genotypes/ukb_la_snps.vcf")
DEFAULT_OUTPUT_DIR = Path("analysis")

app = typer.Typer(
    add_completion=False,
    help=(
        "Activity 2.1.11.1 — join mock UKB phenotype CSV and LA-SNP VCF on eid, "
        "run dominant-model association scans. Synthetic data only."
    ),
)


@app.command()
def integrate(
    pheno: Path = typer.Option(
        DEFAULT_PHENO,
        "--pheno",
        path_type=Path,
        help="Mock UKB phenotype CSV (eid + v2 fields).",
    ),
    vcf: Path = typer.Option(
        DEFAULT_VCF,
        "--vcf",
        path_type=Path,
        help="Mock LA-SNP VCF with sample IDs equal to eid.",
    ),
    output_dir: Path = typer.Option(
        DEFAULT_OUTPUT_DIR,
        "--output-dir",
        path_type=Path,
        help="Directory for assoc_la_snp_*.csv outputs.",
    ),
    verbose: bool = typer.Option(
        False,
        "-v",
        "--verbose",
        help="Enable debug logging.",
    ),
) -> int:
    """Join mock UKB phenotype CSV and LA-SNP VCF; run association scans.

    Args:
        pheno: Mock UKB phenotype CSV (``eid`` + v2 fields).
        vcf: Mock LA-SNP VCF with sample IDs equal to ``eid``.
        output_dir: Directory for ``assoc_la_snp_*.csv`` outputs.
        verbose: When true, enable debug logging.

    Returns:
        Process exit code (``0`` on success, ``1`` if inputs are missing).
    """
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
    if not pheno.is_file():
        print(f"Phenotype file not found: {pheno}", file=sys.stderr)
        return 1
    if not vcf.is_file():
        print(f"VCF not found: {vcf}", file=sys.stderr)
        return 1

    joined, parental, ad = run_integration_pipeline(pheno, vcf, output_dir)
    print(f"Activity {ACTIVITY_ID} — {SYNTHETIC_DISCLAIMER}")
    print(f"Joined cohort: {joined.height} rows")
    print(f"Parental longevity associations: {parental.height} SNPs → {output_dir}")
    print(f"AD diagnosis associations: {ad.height} SNPs → {output_dir}")
    return 0


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for the synthetic UKB integration pipeline.

    Args:
        argv: Optional argument list; ``None`` uses ``sys.argv``.

    Returns:
        Process exit code.
    """
    kwargs: dict[str, Any] = {"standalone_mode": False}
    if argv is not None:
        kwargs["args"] = list(argv)
    try:
        result = app(**kwargs)
    except typer.Exit as exc:
        return int(exc.exit_code)
    except SystemExit as exc:
        code = exc.code
        if code is None:
            return 0
        if isinstance(code, int):
            return code
        return 1
    return 0 if result is None else int(result)
