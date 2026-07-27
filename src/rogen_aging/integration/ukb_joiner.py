"""Deprecated alias — use :mod:`rogen_aging.ukb_integration.ukb_joiner`."""

from __future__ import annotations

import warnings

from rogen_aging.ukb_integration.ukb_joiner import (
    ACTIVITY_ID,
    AD_OUT,
    LA_SNP_ASSOC_COLUMNS,
    PARENTAL_LONGEVITY_OUT,
    SYNTHETIC_DISCLAIMER,
    ad_diagnosis_from_code,
    dominant_odds_ratio,
    genotype_phenotype_contingency,
    join_phenotypes_genotypes,
    load_genotype_matrix_from_vcf,
    load_phenotype_table,
    run_association_scan,
    run_integration_pipeline,
    write_association_results,
)

warnings.warn(
    "rogen_aging.integration.ukb_joiner is deprecated; use rogen_aging.ukb_integration.ukb_joiner",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = [
    "ACTIVITY_ID",
    "AD_OUT",
    "LA_SNP_ASSOC_COLUMNS",
    "PARENTAL_LONGEVITY_OUT",
    "SYNTHETIC_DISCLAIMER",
    "ad_diagnosis_from_code",
    "dominant_odds_ratio",
    "genotype_phenotype_contingency",
    "join_phenotypes_genotypes",
    "load_genotype_matrix_from_vcf",
    "load_phenotype_table",
    "run_association_scan",
    "run_integration_pipeline",
    "write_association_results",
]
