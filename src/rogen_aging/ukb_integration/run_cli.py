"""CLI for Activity 2.1.11.1 synthetic UKB integrative validation."""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any

import typer

from rogen_aging.config import cfg_path, get_config
from rogen_aging.config.cli import config_option, load_cli_config
from rogen_aging.ukb_integration.ukb_joiner import (
    ACTIVITY_ID,
    DEFAULT_AUDIT_LOG,
    MAX_JOIN_DROP_RATE,
    SYNTHETIC_DISCLAIMER,
    JoinDropRateError,
    max_join_drop_rate,
    run_integration_pipeline,
)

app = typer.Typer(
    add_completion=False,
    help=(
        "Activity 2.1.11.1 — join mock UKB phenotype CSV and LA-SNP VCF on eid, "
        "run dominant-model association scans. Synthetic data only."
    ),
)


@app.command()
def integrate(
    config: Path | None = config_option(),
    pheno: Path | None = typer.Option(
        None,
        "--pheno",
        help="Mock UKB phenotype CSV (eid + v2 fields). Default: from config.",
    ),
    vcf: Path | None = typer.Option(
        None,
        "--vcf",
        help="Mock LA-SNP VCF with sample IDs equal to eid. Default: from config.",
    ),
    output_dir: Path | None = typer.Option(
        None,
        "--output-dir",
        help="Directory for assoc_la_snp_*.csv outputs. Default: from config.",
    ),
    audit_log: Path | None = typer.Option(
        None,
        "--audit-log",
        help="Log path for dropped records and eid schema mismatches. Default: from config.",
    ),
    max_drop_rate: float | None = typer.Option(
        None,
        "--max-drop-rate",
        min=0.0,
        max=1.0,
        help="Halt if unmatched eid fraction of the ID union exceeds this value. Default: from config.",
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
        config: Optional YAML config override.
        pheno: Mock UKB phenotype CSV (``eid`` + v2 fields).
        vcf: Mock LA-SNP VCF with sample IDs equal to ``eid``.
        output_dir: Directory for ``assoc_la_snp_*.csv`` outputs.
        audit_log: Destination for dropped-ID / schema mismatch audit logging.
        max_drop_rate: Maximum allowed unmatched ``eid`` fraction before halt.
        verbose: When true, enable debug logging.

    Returns:
        Process exit code (``0`` on success, ``1`` on missing inputs or excess drops).
    """
    load_cli_config(config)
    cfg = get_config()
    resolved_pheno = pheno or cfg_path(cfg, "paths", "ukb", "mock_pheno")
    resolved_vcf = vcf or cfg_path(cfg, "paths", "ukb", "mock_vcf")
    resolved_output = output_dir or cfg_path(cfg, "paths", "ukb", "integration_output_dir")
    resolved_audit = audit_log or cfg_path(cfg, "paths", "ukb", "audit_log")
    resolved_drop = float(max_join_drop_rate() if max_drop_rate is None else max_drop_rate)

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
    if not resolved_pheno.is_file():
        print(f"Phenotype file not found: {resolved_pheno}", file=sys.stderr)
        return 1
    if not resolved_vcf.is_file():
        print(f"VCF not found: {resolved_vcf}", file=sys.stderr)
        return 1

    try:
        joined, parental, ad = run_integration_pipeline(
            resolved_pheno,
            resolved_vcf,
            resolved_output,
            audit_log=resolved_audit,
            max_drop_rate=resolved_drop,
        )
    except JoinDropRateError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Activity {ACTIVITY_ID} — {SYNTHETIC_DISCLAIMER}")
    print(f"Joined cohort: {joined.height} rows")
    print(f"Parental longevity associations: {parental.height} SNPs → {resolved_output}")
    print(f"AD diagnosis associations: {ad.height} SNPs → {resolved_output}")
    print(f"eid join audit log: {resolved_audit}")
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


# Re-export defaults for tests / callers that import module-level paths.
DEFAULT_PHENO = Path("test_data/mock_ukb_rap/phenotypes/ukb_phenotypes.csv")
DEFAULT_VCF = Path("test_data/mock_ukb_rap/genotypes/ukb_la_snps.vcf")
DEFAULT_OUTPUT_DIR = Path("analysis")

__all__ = [
    "DEFAULT_AUDIT_LOG",
    "DEFAULT_OUTPUT_DIR",
    "DEFAULT_PHENO",
    "DEFAULT_VCF",
    "MAX_JOIN_DROP_RATE",
    "app",
    "integrate",
    "main",
]
