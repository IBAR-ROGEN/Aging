"""Deprecated alias for :mod:`rogen_aging.ukb_integration`.

Import ``rogen_aging.ukb_integration`` instead. This package remains as a
compatibility shim so existing ``from rogen_aging.integration …`` imports
continue to work during migration.

Note: ``rogen_aging.ukb_integration`` (synthetic UKB join) is distinct from
``rogen_aging.integrative`` (variant×tissue×phenotype risk scoring).
"""

from __future__ import annotations

import warnings

from rogen_aging.ukb_integration import (
    LA_SNP_ASSOC_COLUMNS,
    ad_diagnosis_from_code,
    join_phenotypes_genotypes,
    load_genotype_matrix_from_vcf,
    load_phenotype_table,
    run_association_scan,
    run_integration_pipeline,
    write_association_results,
)

warnings.warn(
    "rogen_aging.integration is deprecated; use rogen_aging.ukb_integration",
    DeprecationWarning,
    stacklevel=2,
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
