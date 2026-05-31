#!/usr/bin/env python3
"""Activity 2.1.11.1 — thin CLI for synthetic UKB integrative validation.

Joins mock RAP phenotype CSV and LA-SNP VCF (from ``scripts/ukb_mock_gen.py``), runs
dominant-model association scans, and writes tidy result tables. **Synthetic data
only**; no biological conclusions.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from rogen_aging.integration.ukb_joiner import (
    ACTIVITY_ID,
    SYNTHETIC_DISCLAIMER,
    run_integration_pipeline,
)

DEFAULT_PHENO = Path("test_data/mock_ukb_rap/phenotypes/ukb_phenotypes.csv")
DEFAULT_VCF = Path("test_data/mock_ukb_rap/genotypes/ukb_la_snps.vcf")
DEFAULT_OUTPUT_DIR = Path("analysis")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--pheno",
        type=Path,
        default=DEFAULT_PHENO,
        help="Mock UKB phenotype CSV (eid + v2 fields).",
    )
    parser.add_argument(
        "--vcf",
        type=Path,
        default=DEFAULT_VCF,
        help="Mock LA-SNP VCF with sample IDs equal to eid.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for assoc_la_snp_*.csv outputs.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
    if not args.pheno.is_file():
        print(f"Phenotype file not found: {args.pheno}", file=sys.stderr)
        return 1
    if not args.vcf.is_file():
        print(f"VCF not found: {args.vcf}", file=sys.stderr)
        return 1

    joined, parental, ad = run_integration_pipeline(args.pheno, args.vcf, args.output_dir)
    print(f"Activity {ACTIVITY_ID} — {SYNTHETIC_DISCLAIMER}")
    print(f"Joined cohort: {joined.height} rows")
    print(f"Parental longevity associations: {parental.height} SNPs → {args.output_dir}")
    print(f"AD diagnosis associations: {ad.height} SNPs → {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
