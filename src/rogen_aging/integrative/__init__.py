"""Integrative multi-omics analysis: variants × tissues × phenotypic risk.

Consolidates standalone GTEx / VEP / AlphaGenome annotation and phenotype-risk
scripts into a reusable package API. The primary entry point is
:func:`run_integrative_pipeline`.
"""

from __future__ import annotations

from typing import Any

import polars as pl

from rogen_aging.integrative.phenotype_integrator import (
    ALPHAMISSENSE_HIGH_THRESHOLD,
    DEFAULT_WEIGHTS,
    VEP_IMPACT_SCORES,
    PhenotypeIntegrator,
)
from rogen_aging.integrative.variant_tissue_mapper import (
    DEFAULT_TARGET_TISSUES,
    VariantTissueMapper,
    gtex_chromosome,
    normalize_chrom,
    normalize_rsid,
    variant_key,
)


def run_integrative_pipeline(
    variants: pl.DataFrame,
    eqtls: pl.DataFrame,
    *,
    alphagenome: pl.DataFrame | None = None,
    probe_annotation: pl.DataFrame | None = None,
    sample_phenotypes: pl.DataFrame | None = None,
    target_tissues: tuple[str, ...] | list[str] | None = None,
    risk_weights: dict[str, float] | None = None,
    age_accel_col: str = "age_acceleration",
    genotype_col: str = "alt_dosage",
    sample_id_col: str = "sample_id",
    variant_id_col: str = "rsid",
) -> dict[str, Any]:
    """Run the full integrative variant→tissue→phenotype pipeline.

    Cross-references VEP/AlphaGenome-annotated variants with tissue-specific
    GTEx eQTL profiles, optionally links epigenetic methylation probes, and
    computes composite phenotypic risk scores.

    Args:
        variants: Annotated variant table (``chrom``, ``pos``, ``ref``, ``alt``,
            plus optional VEP / score columns).
        eqtls: Long eQTL table (``rsid``, ``tissue``, ``nes``, ``p_value``).
        alphagenome: Optional AlphaGenome score matrix.
        probe_annotation: Optional HM450/EPIC probe→gene annotation.
        sample_phenotypes: Optional long genotype table for sample-level risk.
        target_tissues: GTEx tissues retained by the mapper.
        risk_weights: Optional overrides for :class:`PhenotypeIntegrator` weights.
        age_accel_col: Epigenetic age-acceleration column on variants.
        genotype_col: Dosage column on ``sample_phenotypes``.
        sample_id_col: Sample id column on ``sample_phenotypes``.
        variant_id_col: Shared variant id for sample aggregation.

    Returns:
        Dictionary containing:

        * ``annotated`` — variants with GTEx (+ optional AlphaGenome) columns
        * ``eqtl_summary`` — per-rsID eQTL summary
        * ``variant_risks`` — annotated variants with ``composite_risk``
        * ``methylation_links`` — present when ``probe_annotation`` is given
        * ``sample_profiles`` — present when ``sample_phenotypes`` is given
        * ``mapper`` / ``integrator`` — the configured helper instances
    """
    mapper = VariantTissueMapper(target_tissues=target_tissues)
    mapped = mapper.map_variants_to_tissues(
        variants,
        eqtls,
        alphagenome=alphagenome,
        probe_annotation=probe_annotation,
    )

    integrator = PhenotypeIntegrator(weights=risk_weights)
    risk = integrator.build_risk_profile(
        mapped["annotated"],
        sample_phenotypes=sample_phenotypes,
        age_accel_col=age_accel_col,
        genotype_col=genotype_col,
        sample_id_col=sample_id_col,
        variant_id_col=variant_id_col,
    )

    result: dict[str, Any] = {
        "annotated": mapped["annotated"],
        "eqtl_summary": mapped["eqtl_summary"],
        "variant_risks": risk["variant_risks"],
        "mapper": mapper,
        "integrator": integrator,
    }
    if "methylation_links" in mapped:
        result["methylation_links"] = mapped["methylation_links"]
    if "sample_profiles" in risk:
        result["sample_profiles"] = risk["sample_profiles"]
    return result


__all__ = [
    "ALPHAMISSENSE_HIGH_THRESHOLD",
    "DEFAULT_TARGET_TISSUES",
    "DEFAULT_WEIGHTS",
    "PhenotypeIntegrator",
    "VEP_IMPACT_SCORES",
    "VariantTissueMapper",
    "gtex_chromosome",
    "normalize_chrom",
    "normalize_rsid",
    "run_integrative_pipeline",
    "variant_key",
]
