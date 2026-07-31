"""Link molecular variant scores with composite phenotypic risk profiles.

Combines VEP impact, AlphaGenome / AlphaMissense scores, GTEx eQTL evidence,
and optional epigenetic (methylation age-acceleration) markers into a single
per-variant or per-sample risk table.
"""

from __future__ import annotations

import polars as pl

from rogen_aging.config import (
    alphamissense_high_threshold,
    risk_weights,
    vep_impact_scores,
)

# Module-level aliases kept for public API / test imports. Values are seeded from
# config/default.yaml via :func:`rogen_aging.config.get_config` and refreshed when
# :func:`rogen_aging.config.set_config` / ``load_config`` runs.
DEFAULT_WEIGHTS: dict[str, float] = risk_weights()
VEP_IMPACT_SCORES: dict[str, float] = vep_impact_scores()
ALPHAMISSENSE_HIGH_THRESHOLD: float = alphamissense_high_threshold()


class PhenotypeIntegrator:
    """Integrate molecular annotation scores into composite phenotypic risk.

    Attributes:
        weights: Relative contribution of each evidence channel. Values are
            renormalised to sum to 1.0 at score time.
    """

    def __init__(self, weights: dict[str, float] | None = None) -> None:
        """Initialise the integrator.

        Args:
            weights: Optional overrides for channel weights. Missing keys fall
                back to configured :data:`DEFAULT_WEIGHTS`.
        """
        merged = dict(risk_weights())
        if weights is not None:
            merged.update(weights)
        self.weights: dict[str, float] = merged

    def score_vep_impact(self, frame: pl.DataFrame) -> pl.Series:
        """Map VEP impact labels onto a ``[0, 1]`` severity score.

        Args:
            frame: Table with optional ``vep_impact`` column.

        Returns:
            Float series aligned to ``frame`` rows (0.0 when impact is missing).
        """
        if "vep_impact" not in frame.columns:
            return pl.Series("vep_impact_score", [0.0] * frame.height)
        scores = vep_impact_scores()
        mapping = pl.DataFrame(
            {
                "vep_impact": list(scores.keys()),
                "vep_impact_score": list(scores.values()),
            }
        )
        joined = frame.select(
            pl.col("vep_impact").cast(pl.Utf8).str.to_uppercase().alias("vep_impact")
        ).join(mapping, on="vep_impact", how="left")
        return joined["vep_impact_score"].fill_null(0.0)

    def score_alphagenome(self, frame: pl.DataFrame) -> pl.Series:
        """Score regulatory impact from Absolute AlphaGenome percent change.

        Uses ``alphagenome_abs_perc_change`` when present, otherwise the absolute
        value of ``alphagenome_diff`` / ``alphagenome_perc_change``. Values are
        clipped into ``[0, 1]`` via ``min(abs(x) / 50, 1)`` for percent-change
        inputs and ``min(abs(x), 1)`` for raw diffs.

        Args:
            frame: Table with optional AlphaGenome score columns.

        Returns:
            Float series in ``[0, 1]``.
        """
        n = frame.height
        if "alphagenome_abs_perc_change" in frame.columns:
            raw = frame["alphagenome_abs_perc_change"].cast(pl.Float64, strict=False)
            return (raw.fill_null(0.0).abs() / 50.0).clip(0.0, 1.0).alias("alphagenome_score")
        if "alphagenome_perc_change" in frame.columns:
            raw = frame["alphagenome_perc_change"].cast(pl.Float64, strict=False)
            return (raw.fill_null(0.0).abs() / 50.0).clip(0.0, 1.0).alias("alphagenome_score")
        if "alphagenome_diff" in frame.columns:
            raw = frame["alphagenome_diff"].cast(pl.Float64, strict=False)
            return raw.fill_null(0.0).abs().clip(0.0, 1.0).alias("alphagenome_score")
        return pl.Series("alphagenome_score", [0.0] * n)

    def score_alphamissense(self, frame: pl.DataFrame) -> pl.Series:
        """Return AlphaMissense pathogenicity scores clipped to ``[0, 1]``.

        Args:
            frame: Table with optional ``alphamissense_score`` column.

        Returns:
            Float series in ``[0, 1]`` (0.0 when missing).
        """
        if "alphamissense_score" not in frame.columns:
            return pl.Series("alphamissense_channel", [0.0] * frame.height)
        raw = frame["alphamissense_score"].cast(pl.Float64, strict=False)
        return raw.fill_null(0.0).clip(0.0, 1.0).alias("alphamissense_channel")

    def score_gtex_eqtl(self, frame: pl.DataFrame) -> pl.Series:
        """Score tissue eQTL evidence from hit count and best p-value.

        Combines a saturating hit-count term with ``-log10(p)`` of the best
        eQTL, each mapped into ``[0, 1]`` and averaged when both exist.

        Args:
            frame: Table with optional ``gtex_n_eqtls`` and ``gtex_best_p_value``.

        Returns:
            Float series in ``[0, 1]``.
        """
        n = frame.height
        hit_term = pl.Series("hit", [0.0] * n)
        if "gtex_n_eqtls" in frame.columns:
            counts = frame["gtex_n_eqtls"].cast(pl.Float64, strict=False).fill_null(0.0)
            hit_term = (counts / (counts + 5.0)).clip(0.0, 1.0)

        p_term = pl.Series("p", [0.0] * n)
        if "gtex_best_p_value" in frame.columns:
            pvals = frame["gtex_best_p_value"].cast(pl.Float64, strict=False)
            # Avoid log of zero/negative; treat null as no evidence.
            safe = pvals.fill_null(1.0).clip(1e-300, 1.0)
            logp = (-safe.log10() / 10.0).clip(0.0, 1.0)
            p_term = logp

        if "gtex_n_eqtls" in frame.columns and "gtex_best_p_value" in frame.columns:
            return ((hit_term + p_term) / 2.0).alias("gtex_eqtl_score")
        if "gtex_n_eqtls" in frame.columns:
            return hit_term.alias("gtex_eqtl_score")
        if "gtex_best_p_value" in frame.columns:
            return p_term.alias("gtex_eqtl_score")
        return pl.Series("gtex_eqtl_score", [0.0] * n)

    def score_epigenetic(
        self,
        frame: pl.DataFrame,
        *,
        age_accel_col: str = "age_acceleration",
    ) -> pl.Series:
        """Score epigenetic burden from methylation age acceleration.

        Absolute age acceleration (years) is mapped via
        ``min(|accel| / 10, 1)``. When the column is absent, returns zeros.

        Args:
            frame: Table that may contain age-acceleration values.
            age_accel_col: Column name for DNAm age − chronological age.

        Returns:
            Float series in ``[0, 1]``.
        """
        if age_accel_col not in frame.columns:
            return pl.Series("epigenetic_score", [0.0] * frame.height)
        raw = frame[age_accel_col].cast(pl.Float64, strict=False).fill_null(0.0)
        return (raw.abs() / 10.0).clip(0.0, 1.0).alias("epigenetic_score")

    def compute_composite_risk(
        self,
        annotated_variants: pl.DataFrame,
        *,
        age_accel_col: str = "age_acceleration",
    ) -> pl.DataFrame:
        """Attach channel scores and a weighted composite risk column.

        Args:
            annotated_variants: Variants already joined with VEP, AlphaGenome,
                AlphaMissense, and/or GTEx summary columns.
            age_accel_col: Optional epigenetic age-acceleration column name.

        Returns:
            Input columns plus ``*_score`` channels and ``composite_risk``.
        """
        frame = annotated_variants
        vep = self.score_vep_impact(frame)
        ag = self.score_alphagenome(frame)
        am = self.score_alphamissense(frame)
        gtex = self.score_gtex_eqtl(frame)
        epi = self.score_epigenetic(frame, age_accel_col=age_accel_col)

        scored = frame.with_columns(
            vep.alias("vep_impact_score"),
            ag.alias("alphagenome_score"),
            am.alias("alphamissense_channel"),
            gtex.alias("gtex_eqtl_score"),
            epi.alias("epigenetic_score"),
        )

        total_w = sum(self.weights.values())
        if total_w <= 0:
            raise ValueError("PhenotypeIntegrator weights must sum to a positive value")

        composite = (
            scored["vep_impact_score"] * (self.weights["vep_impact"] / total_w)
            + scored["alphagenome_score"] * (self.weights["alphagenome"] / total_w)
            + scored["alphamissense_channel"] * (self.weights["alphamissense"] / total_w)
            + scored["gtex_eqtl_score"] * (self.weights["gtex_eqtl"] / total_w)
            + scored["epigenetic_score"] * (self.weights["epigenetic"] / total_w)
        )
        return scored.with_columns(composite.alias("composite_risk"))

    def integrate_sample_profiles(
        self,
        variant_risks: pl.DataFrame,
        sample_phenotypes: pl.DataFrame,
        *,
        genotype_col: str = "alt_dosage",
        sample_id_col: str = "sample_id",
        variant_id_col: str = "rsid",
    ) -> pl.DataFrame:
        """Aggregate per-variant composite risk into per-sample risk profiles.

        Expects a long genotype table (``sample_id``, ``rsid``, dosage) joined
        against per-variant ``composite_risk`` scores. Sample-level risk is the
        dosage-weighted mean of variant composite scores, plus any phenotypic
        columns carried on ``sample_phenotypes``.

        Args:
            variant_risks: Output of :meth:`compute_composite_risk` (needs
                ``composite_risk`` and ``variant_id_col``).
            sample_phenotypes: Long genotypes with optional phenotype columns
                already joined, or a sample metadata table joined after
                aggregation when genotypes are supplied separately via the same
                frame containing ``genotype_col``.
            genotype_col: Allele dosage column (0/1/2).
            sample_id_col: Sample identifier column.
            variant_id_col: Variant identifier shared with ``variant_risks``.

        Returns:
            One row per sample with ``sample_risk`` and phenotype columns.

        Raises:
            ValueError: If required columns are missing or ``composite_risk``
                has not been computed.
        """
        if "composite_risk" not in variant_risks.columns:
            raise ValueError(
                "variant_risks must include composite_risk; call compute_composite_risk first"
            )
        for col, label in (
            (sample_id_col, "sample_phenotypes"),
            (variant_id_col, "sample_phenotypes"),
            (genotype_col, "sample_phenotypes"),
        ):
            if col not in sample_phenotypes.columns:
                raise ValueError(f"{label} missing required column: {col}")
        if variant_id_col not in variant_risks.columns:
            raise ValueError(f"variant_risks missing column: {variant_id_col}")

        risk_cols = variant_risks.select([variant_id_col, "composite_risk"]).unique(
            subset=[variant_id_col], keep="first"
        )

        joined = sample_phenotypes.join(risk_cols, on=variant_id_col, how="left")
        joined = joined.with_columns(
            (
                pl.col(genotype_col).cast(pl.Float64, strict=False).fill_null(0.0)
                * pl.col("composite_risk").fill_null(0.0)
            ).alias("_weighted_risk"),
            pl.col(genotype_col).cast(pl.Float64, strict=False).fill_null(0.0).alias("_dosage"),
        )

        phenotype_cols = [
            c
            for c in sample_phenotypes.columns
            if c not in {sample_id_col, variant_id_col, genotype_col}
        ]

        agg_exprs: list[pl.Expr] = [
            pl.col("_weighted_risk").sum().alias("_risk_sum"),
            pl.col("_dosage").sum().alias("_dosage_sum"),
            pl.len().alias("n_variants"),
        ]
        for col in phenotype_cols:
            agg_exprs.append(pl.col(col).first().alias(col))

        per_sample = joined.group_by(sample_id_col).agg(agg_exprs)
        return per_sample.with_columns(
            pl.when(pl.col("_dosage_sum") > 0)
            .then(pl.col("_risk_sum") / pl.col("_dosage_sum"))
            .otherwise(0.0)
            .alias("sample_risk")
        ).drop(["_risk_sum", "_dosage_sum"])

    def build_risk_profile(
        self,
        annotated_variants: pl.DataFrame,
        sample_phenotypes: pl.DataFrame | None = None,
        *,
        age_accel_col: str = "age_acceleration",
        genotype_col: str = "alt_dosage",
        sample_id_col: str = "sample_id",
        variant_id_col: str = "rsid",
    ) -> dict[str, pl.DataFrame]:
        """Compute variant-level and optional sample-level risk profiles.

        Args:
            annotated_variants: Molecularly annotated variants.
            sample_phenotypes: Optional long genotype (+ phenotype) table.
            age_accel_col: Epigenetic age-acceleration column on variants.
            genotype_col: Dosage column on ``sample_phenotypes``.
            sample_id_col: Sample id column on ``sample_phenotypes``.
            variant_id_col: Shared variant id column.

        Returns:
            Dictionary with ``variant_risks`` and, when samples are provided,
            ``sample_profiles``.
        """
        variant_risks = self.compute_composite_risk(annotated_variants, age_accel_col=age_accel_col)
        result: dict[str, pl.DataFrame] = {"variant_risks": variant_risks}
        if sample_phenotypes is not None:
            result["sample_profiles"] = self.integrate_sample_profiles(
                variant_risks,
                sample_phenotypes,
                genotype_col=genotype_col,
                sample_id_col=sample_id_col,
                variant_id_col=variant_id_col,
            )
        return result


__all__ = [
    "ALPHAMISSENSE_HIGH_THRESHOLD",
    "DEFAULT_WEIGHTS",
    "PhenotypeIntegrator",
    "VEP_IMPACT_SCORES",
]
