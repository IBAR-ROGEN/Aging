"""Synthetic UKB integrative validation (Activity 2.1.11.1); not for real cohort data."""

from __future__ import annotations

from rogen_aging.integration.ukb_joiner import (
    LA_SNP_ASSOC_COLUMNS,
    ad_diagnosis_from_code,
    join_phenotypes_genotypes,
    load_genotype_matrix_from_vcf,
    load_phenotype_table,
    run_association_scan,
    run_integration_pipeline,
    write_association_results,
)

__all__ = [
    "LA_SNP_ASSOC_COLUMNS",
    "ad_diagnosis_from_code",
    "join_phenotypes_genotypes",
    "load_genotype_matrix_from_vcf",
    "load_phenotype_table",
    "run_association_scan",
    "run_integration_pipeline",
    "write_association_results",
]
