"""Cross-reference VEP/AlphaGenome annotated variants with tissue eQTL profiles.

Consolidates join and summarisation logic previously scattered across standalone
GTEx / VEP annotation scripts (``scripts/ukb/annotate_la_snps_gtex.py``,
``run_july_annotation_pipeline.py``) into a reusable, API-free mapper.
"""

from __future__ import annotations

import re
from typing import Final

import polars as pl

DEFAULT_TARGET_TISSUES: Final[tuple[str, ...]] = (
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra",
    "Whole_Blood",
)

RSID_RE: Final[re.Pattern[str]] = re.compile(r"^rs\d+$", re.IGNORECASE)

REQUIRED_VARIANT_COLS: Final[tuple[str, ...]] = (
    "chrom",
    "pos",
    "ref",
    "alt",
)

EQTL_REQUIRED_COLS: Final[tuple[str, ...]] = (
    "rsid",
    "tissue",
    "nes",
    "p_value",
)

METHYLATION_PROBE_COLS: Final[tuple[str, ...]] = (
    "IlmnID",
    "UCSC_RefGene_Name",
)


def normalize_chrom(value: object) -> str:
    """Return chromosome without a ``chr`` prefix.

    Args:
        value: Raw chromosome label (e.g. ``\"chr17\"`` or ``17``).

    Returns:
        Chromosome string without a leading ``chr`` prefix.
    """
    text = str(value).strip()
    if text.lower().startswith("chr"):
        return text[3:]
    return text


def gtex_chromosome(value: object) -> str:
    """Format chromosome as GTEx expects (``chrN``).

    Args:
        value: Raw chromosome label.

    Returns:
        Chromosome string with a ``chr`` prefix.
    """
    return f"chr{normalize_chrom(value)}"


def normalize_rsid(raw: object) -> str | None:
    """Return a normalized rsID, or ``None`` when the value is not rs-formatted.

    Args:
        raw: Candidate rsID value.

    Returns:
        Lowercase-preserving rsID string, or ``None`` if invalid.
    """
    if raw is None:
        return None
    token = str(raw).strip()
    if not token or token.lower() in {"nan", "none", "null"}:
        return None
    if not RSID_RE.match(token):
        return None
    return token


def variant_key(chrom: object, pos: object, ref: object, alt: object) -> str:
    """Build a locus key used for joining annotation and score matrices.

    Args:
        chrom: Chromosome.
        pos: Genomic position (1-based).
        ref: Reference allele.
        alt: Alternate allele (first allele used if multi-allelic).

    Returns:
        Key of the form ``chrom:pos:REF:ALT``.
    """
    primary_alt = str(alt).split(",")[0].strip()
    position = int(str(pos))
    return (
        f"{normalize_chrom(chrom)}:{position}:" f"{str(ref).strip().upper()}:{primary_alt.upper()}"
    )


class VariantTissueMapper:
    """Map VEP/AlphaGenome-annotated variants onto tissue-specific eQTL profiles.

    The mapper operates on in-memory Polars frames so unit tests and offline
    pipelines never require live GTEx or Ensembl API calls. Optional methylation
    probe annotations can be linked via shared gene symbols.

    Attributes:
        target_tissues: GTEx ``tissueSiteDetailId`` values retained in summaries.
    """

    def __init__(
        self,
        target_tissues: tuple[str, ...] | list[str] | None = None,
    ) -> None:
        """Initialise the mapper.

        Args:
            target_tissues: Tissues to retain when filtering eQTL tables.
                Defaults to brain regions plus whole blood.
        """
        tissues = target_tissues if target_tissues is not None else DEFAULT_TARGET_TISSUES
        self.target_tissues: tuple[str, ...] = tuple(tissues)
        self._target_set: frozenset[str] = frozenset(self.target_tissues)

    def prepare_variants(self, variants: pl.DataFrame) -> pl.DataFrame:
        """Normalise variant coordinates and attach a ``variant_key`` column.

        Args:
            variants: Frame with at least ``chrom``, ``pos``, ``ref``, ``alt``.
                Optional columns ``rsid`` / ``gene_symbol`` are preserved.

        Returns:
            Copy with normalised alleles, optional ``rsid``, and ``variant_key``.

        Raises:
            ValueError: If required coordinate columns are missing.
        """
        missing = [c for c in REQUIRED_VARIANT_COLS if c not in variants.columns]
        if missing:
            raise ValueError(f"Variants frame missing required columns: {missing}")

        frame = variants.with_columns(
            pl.col("chrom").cast(pl.Utf8).map_elements(normalize_chrom, return_dtype=pl.Utf8),
            pl.col("pos").cast(pl.Int64),
            pl.col("ref").cast(pl.Utf8).str.strip_chars().str.to_uppercase(),
            pl.col("alt")
            .cast(pl.Utf8)
            .str.split(",")
            .list.first()
            .str.strip_chars()
            .str.to_uppercase()
            .alias("alt"),
        )
        if "rsid" in frame.columns:
            frame = frame.with_columns(
                pl.col("rsid").map_elements(normalize_rsid, return_dtype=pl.Utf8).alias("rsid")
            )
        frame = frame.with_columns(
            (
                pl.col("chrom")
                + pl.lit(":")
                + pl.col("pos").cast(pl.Utf8)
                + pl.lit(":")
                + pl.col("ref")
                + pl.lit(":")
                + pl.col("alt")
            ).alias("variant_key")
        )
        return frame

    def filter_eqtls_to_target_tissues(self, eqtls: pl.DataFrame) -> pl.DataFrame:
        """Keep eQTL rows whose tissue is in ``target_tissues``.

        Args:
            eqtls: Long eQTL table with a ``tissue`` column.

        Returns:
            Filtered frame (may be empty).

        Raises:
            ValueError: If ``tissue`` is missing.
        """
        if "tissue" not in eqtls.columns:
            raise ValueError("eQTL frame missing required column: tissue")
        return eqtls.filter(pl.col("tissue").is_in(list(self._target_set)))

    def summarize_eqtl_profiles(self, eqtls: pl.DataFrame) -> pl.DataFrame:
        """Collapse long eQTL hits into one summary row per variant rsID.

        Args:
            eqtls: Long table with ``rsid``, ``tissue``, ``nes``, ``p_value``.
                Optional ``gene_symbol`` / ``eqtl_gene_symbol`` and
                ``gtex_variant_id`` are summarised when present.

        Returns:
            Wide summary with best (lowest p) tissue hit and tissue list.

        Raises:
            ValueError: If required eQTL columns are missing.
        """
        missing = [c for c in EQTL_REQUIRED_COLS if c not in eqtls.columns]
        if missing:
            raise ValueError(f"eQTL frame missing required columns: {missing}")

        filtered = self.filter_eqtls_to_target_tissues(eqtls)
        if filtered.is_empty():
            return pl.DataFrame(
                schema={
                    "rsid": pl.Utf8,
                    "gtex_n_eqtls": pl.Int64,
                    "gtex_best_tissue": pl.Utf8,
                    "gtex_best_gene": pl.Utf8,
                    "gtex_best_slope": pl.Float64,
                    "gtex_best_p_value": pl.Float64,
                    "gtex_tissues": pl.Utf8,
                    "gtex_variant_id": pl.Utf8,
                }
            )

        gene_col = (
            "eqtl_gene_symbol"
            if "eqtl_gene_symbol" in filtered.columns
            else "gene_symbol" if "gene_symbol" in filtered.columns else None
        )

        ranked = filtered.sort(["rsid", "p_value"], nulls_last=True)
        agg_exprs: list[pl.Expr] = [
            pl.len().alias("gtex_n_eqtls"),
            pl.col("tissue").first().alias("gtex_best_tissue"),
            pl.col("nes").first().alias("gtex_best_slope"),
            pl.col("p_value").first().alias("gtex_best_p_value"),
            pl.col("tissue").unique().sort().str.join(";").alias("gtex_tissues"),
        ]
        if gene_col is not None:
            agg_exprs.insert(2, pl.col(gene_col).first().alias("gtex_best_gene"))
        best = ranked.group_by("rsid").agg(agg_exprs)
        if gene_col is None:
            best = best.with_columns(pl.lit(None).cast(pl.Utf8).alias("gtex_best_gene"))

        if "gtex_variant_id" in filtered.columns:
            ids = (
                filtered.select(["rsid", "gtex_variant_id"])
                .drop_nulls(subset=["gtex_variant_id"])
                .unique(subset=["rsid"], keep="first")
            )
            best = best.join(ids, on="rsid", how="left")
        else:
            best = best.with_columns(pl.lit(None).cast(pl.Utf8).alias("gtex_variant_id"))

        return best

    def join_annotations_with_eqtls(
        self,
        variants: pl.DataFrame,
        eqtls: pl.DataFrame,
        *,
        on: str = "rsid",
    ) -> pl.DataFrame:
        """Left-join annotated variants with per-variant eQTL summaries.

        Args:
            variants: VEP/AlphaGenome annotated variants (will be normalised).
            eqtls: Long eQTL table to summarise before joining.
            on: Join key (``\"rsid\"`` or ``\"variant_key\"``).

        Returns:
            Annotated variants with GTEx summary columns attached.

        Raises:
            ValueError: If the join key is absent after preparation.
        """
        prepared = self.prepare_variants(variants)
        if on not in prepared.columns:
            raise ValueError(f"Prepared variants missing join key: {on}")

        if on == "rsid":
            summary = self.summarize_eqtl_profiles(eqtls)
            return prepared.join(summary, on="rsid", how="left")

        eqtl_frame = eqtls
        if "variant_key" not in eqtl_frame.columns:
            if not all(c in eqtl_frame.columns for c in REQUIRED_VARIANT_COLS):
                raise ValueError(
                    "Cannot join on variant_key: eQTL table lacks coordinates or variant_key"
                )
            eqtl_frame = self.prepare_variants(eqtl_frame)

        select_cols = ["variant_key", "tissue", "nes", "p_value"]
        for optional in ("gene_symbol", "eqtl_gene_symbol", "gtex_variant_id"):
            if optional in eqtl_frame.columns:
                select_cols.append(optional)
        keyed = eqtl_frame.select(select_cols).with_columns(pl.col("variant_key").alias("rsid"))
        key_summary = self.summarize_eqtl_profiles(keyed).rename({"rsid": "variant_key"})
        return prepared.join(key_summary, on="variant_key", how="left")

    def join_alphagenome_scores(
        self,
        variants: pl.DataFrame,
        alphagenome: pl.DataFrame,
    ) -> pl.DataFrame:
        """Left-join AlphaGenome score columns onto annotated variants.

        Prefers ``variant_key`` when both frames expose coordinates; otherwise
        joins on ``rsid``.

        Args:
            variants: Annotated variant table.
            alphagenome: Score matrix with ``alphagenome_*`` columns (aliases
                such as ``ref_score`` / ``diff`` are normalised).

        Returns:
            Variants with AlphaGenome score columns attached.
        """
        prepared = self.prepare_variants(variants)
        scores = self._normalize_alphagenome(alphagenome)

        if "variant_key" in scores.columns:
            keep = [c for c in scores.columns if c.startswith("alphagenome_") or c == "variant_key"]
            return prepared.join(
                scores.select(keep).unique(subset=["variant_key"], keep="first"),
                on="variant_key",
                how="left",
            )
        if "rsid" not in scores.columns or "rsid" not in prepared.columns:
            raise ValueError("AlphaGenome scores need variant_key or rsid to join onto variants")
        keep = [c for c in scores.columns if c.startswith("alphagenome_") or c == "rsid"]
        return prepared.join(
            scores.select(keep).unique(subset=["rsid"], keep="first"),
            on="rsid",
            how="left",
        )

    def map_methylation_markers(
        self,
        variants: pl.DataFrame,
        probe_annotation: pl.DataFrame,
        *,
        gene_col: str = "gene_symbol",
    ) -> pl.DataFrame:
        """Link variants to HM450 / EPIC methylation probes via gene symbols.

        Args:
            variants: Annotated variants with a gene symbol column.
            probe_annotation: Probe table with ``IlmnID`` and
                ``UCSC_RefGene_Name`` (semicolon-separated gene lists allowed).
            gene_col: Gene symbol column on ``variants``.

        Returns:
            Long table of ``variant_key`` / ``rsid`` / gene / ``IlmnID`` links.

        Raises:
            ValueError: If required columns are missing.
        """
        prepared = self.prepare_variants(variants)
        if gene_col not in prepared.columns:
            raise ValueError(f"Variants missing gene column: {gene_col}")
        missing = [c for c in METHYLATION_PROBE_COLS if c not in probe_annotation.columns]
        if missing:
            raise ValueError(f"Probe annotation missing columns: {missing}")

        probes = (
            probe_annotation.select(
                pl.col("IlmnID").cast(pl.Utf8),
                pl.col("UCSC_RefGene_Name").cast(pl.Utf8).alias("gene_list"),
            )
            .filter(pl.col("gene_list").is_not_null() & (pl.col("gene_list") != ""))
            .with_columns(
                pl.col("gene_list").str.split(";").alias("genes"),
            )
            .explode("genes")
            .with_columns(pl.col("genes").str.strip_chars().alias("gene_symbol"))
            .filter(pl.col("gene_symbol") != "")
            .select(["IlmnID", "gene_symbol"])
            .unique()
        )

        var_genes = prepared.select(
            [c for c in ("variant_key", "rsid", gene_col) if c in prepared.columns]
        ).rename({gene_col: "gene_symbol"} if gene_col != "gene_symbol" else {})

        return var_genes.join(probes, on="gene_symbol", how="inner")

    def map_variants_to_tissues(
        self,
        variants: pl.DataFrame,
        eqtls: pl.DataFrame,
        alphagenome: pl.DataFrame | None = None,
        probe_annotation: pl.DataFrame | None = None,
    ) -> dict[str, pl.DataFrame]:
        """Run the full offline variant→tissue (+ optional methylation) map.

        Args:
            variants: VEP/AlphaGenome annotated variants.
            eqtls: Long GTEx eQTL hit table.
            alphagenome: Optional AlphaGenome score matrix.
            probe_annotation: Optional HM450/EPIC probe→gene annotation.

        Returns:
            Dictionary with keys ``annotated``, ``eqtl_summary``, and optionally
            ``methylation_links``.
        """
        annotated = self.join_annotations_with_eqtls(variants, eqtls)
        if alphagenome is not None:
            # Avoid re-normalising twice: join scores onto already prepared frame.
            scores = self._normalize_alphagenome(alphagenome)
            if "variant_key" in scores.columns:
                keep = [
                    c for c in scores.columns if c.startswith("alphagenome_") or c == "variant_key"
                ]
                annotated = annotated.join(
                    scores.select(keep).unique(subset=["variant_key"], keep="first"),
                    on="variant_key",
                    how="left",
                )
            elif "rsid" in scores.columns and "rsid" in annotated.columns:
                keep = [c for c in scores.columns if c.startswith("alphagenome_") or c == "rsid"]
                annotated = annotated.join(
                    scores.select(keep).unique(subset=["rsid"], keep="first"),
                    on="rsid",
                    how="left",
                )

        result: dict[str, pl.DataFrame] = {
            "annotated": annotated,
            "eqtl_summary": self.summarize_eqtl_profiles(eqtls),
        }
        if probe_annotation is not None:
            result["methylation_links"] = self.map_methylation_markers(annotated, probe_annotation)
        return result

    def _normalize_alphagenome(self, frame: pl.DataFrame) -> pl.DataFrame:
        """Rename common AlphaGenome aliases onto canonical column names.

        Args:
            frame: Raw AlphaGenome score matrix.

        Returns:
            Frame with ``alphagenome_*`` columns and join keys when available.
        """
        mapping = {
            "ref_score": "alphagenome_ref_score",
            "alt_score": "alphagenome_alt_score",
            "diff": "alphagenome_diff",
            "perc_change": "alphagenome_perc_change",
            "abs_perc_change": "alphagenome_abs_perc_change",
            "snp": "rsid",
            "gene": "gene_symbol",
        }
        renames = {
            src: dst
            for src, dst in mapping.items()
            if src in frame.columns and dst not in frame.columns
        }
        scores = frame.rename(renames) if renames else frame
        if "rsid" in scores.columns:
            scores = scores.with_columns(
                pl.col("rsid").map_elements(normalize_rsid, return_dtype=pl.Utf8)
            )
        if all(c in scores.columns for c in REQUIRED_VARIANT_COLS):
            scores = self.prepare_variants(scores)
        return scores


__all__ = [
    "DEFAULT_TARGET_TISSUES",
    "VariantTissueMapper",
    "gtex_chromosome",
    "normalize_chrom",
    "normalize_rsid",
    "variant_key",
]
