"""Unit tests for ``rogen_aging.integrative`` (variant×tissue×phenotype)."""

from __future__ import annotations

import polars as pl
import pytest

from rogen_aging.integrative import (
    DEFAULT_TARGET_TISSUES,
    PhenotypeIntegrator,
    VariantTissueMapper,
    gtex_chromosome,
    normalize_chrom,
    normalize_rsid,
    run_integrative_pipeline,
    variant_key,
)
from rogen_aging.integrative.phenotype_integrator import DEFAULT_WEIGHTS, VEP_IMPACT_SCORES
from rogen_aging.integrative.variant_tissue_mapper import (
    normalize_chrom as normalize_chrom_mod,
)


@pytest.fixture
def annotated_variants() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "chrom": ["chr17", "2", "19"],
            "pos": [100, 200, 300],
            "ref": ["A", "G", "C"],
            "alt": ["G", "T,C", "T"],
            "rsid": ["rs1", "rs2", "rs3"],
            "gene_symbol": ["FOXO3", "APOE", "CETP"],
            "vep_impact": ["HIGH", "MODERATE", "LOW"],
            "alphagenome_abs_perc_change": [25.0, 10.0, 0.0],
            "alphamissense_score": [0.8, 0.2, None],
            "age_acceleration": [5.0, -2.0, 0.0],
        }
    )


@pytest.fixture
def eqtl_table() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "rsid": ["rs1", "rs1", "rs2", "rs2", "rs9"],
            "tissue": [
                "Brain_Cortex",
                "Whole_Blood",
                "Brain_Hippocampus",
                "Lung",  # non-target — filtered out
                "Whole_Blood",
            ],
            "nes": [0.5, -0.2, 0.1, 0.9, 0.3],
            "p_value": [1e-8, 1e-4, 1e-3, 1e-10, 1e-5],
            "gene_symbol": ["FOXO3", "FOXO3", "APOE", "APOE", "OTHER"],
            "gtex_variant_id": ["v1", "v1", "v2", "v2", "v9"],
        }
    )


@pytest.fixture
def probe_annotation() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "IlmnID": ["cg0001", "cg0002", "cg0003"],
            "UCSC_RefGene_Name": ["FOXO3", "APOE;TOMM40", "NONE"],
        }
    )


def test_normalize_helpers() -> None:
    assert normalize_chrom("chr17") == "17"
    assert normalize_chrom_mod("17") == "17"
    assert gtex_chromosome("17") == "chr17"
    assert normalize_rsid("rs429358") == "rs429358"
    assert normalize_rsid("not-an-rsid") is None
    assert normalize_rsid(None) is None
    assert variant_key("chr2", 100, "a", "t,g") == "2:100:A:T"


def test_prepare_variants_adds_key(annotated_variants: pl.DataFrame) -> None:
    mapper = VariantTissueMapper()
    prepared = mapper.prepare_variants(annotated_variants)
    assert "variant_key" in prepared.columns
    assert prepared["chrom"].to_list() == ["17", "2", "19"]
    assert prepared["alt"].to_list() == ["G", "T", "T"]
    assert prepared["variant_key"].to_list()[1] == "2:200:G:T"


def test_prepare_variants_missing_cols_raises() -> None:
    mapper = VariantTissueMapper()
    with pytest.raises(ValueError, match="missing required columns"):
        mapper.prepare_variants(pl.DataFrame({"chrom": ["1"]}))


def test_filter_and_summarize_eqtls(eqtl_table: pl.DataFrame) -> None:
    mapper = VariantTissueMapper()
    filtered = mapper.filter_eqtls_to_target_tissues(eqtl_table)
    assert filtered.height == 4  # Lung row dropped
    assert "Lung" not in filtered["tissue"].to_list()

    summary = mapper.summarize_eqtl_profiles(eqtl_table)
    assert set(summary["rsid"].to_list()) == {"rs1", "rs2", "rs9"}
    rs1 = summary.filter(pl.col("rsid") == "rs1").row(0, named=True)
    assert rs1["gtex_n_eqtls"] == 2
    assert rs1["gtex_best_tissue"] == "Brain_Cortex"
    assert rs1["gtex_best_p_value"] == pytest.approx(1e-8)
    assert "Brain_Cortex" in rs1["gtex_tissues"]
    assert "Whole_Blood" in rs1["gtex_tissues"]


def test_join_annotations_with_eqtls(
    annotated_variants: pl.DataFrame, eqtl_table: pl.DataFrame
) -> None:
    mapper = VariantTissueMapper()
    joined = mapper.join_annotations_with_eqtls(annotated_variants, eqtl_table)
    assert joined.height == 3
    rs1 = joined.filter(pl.col("rsid") == "rs1").row(0, named=True)
    assert rs1["gtex_n_eqtls"] == 2
    rs3 = joined.filter(pl.col("rsid") == "rs3").row(0, named=True)
    assert rs3["gtex_n_eqtls"] is None


def test_join_alphagenome_scores(annotated_variants: pl.DataFrame) -> None:
    mapper = VariantTissueMapper()
    scores = pl.DataFrame(
        {
            "snp": ["rs1", "rs2"],
            "diff": [0.4, -0.1],
            "abs_perc_change": [20.0, 5.0],
        }
    )
    joined = mapper.join_alphagenome_scores(
        annotated_variants.select(["chrom", "pos", "ref", "alt", "rsid"]),
        scores,
    )
    assert "alphagenome_diff" in joined.columns
    assert "alphagenome_abs_perc_change" in joined.columns
    assert joined.filter(pl.col("rsid") == "rs1")["alphagenome_diff"][0] == pytest.approx(0.4)


def test_map_methylation_markers(
    annotated_variants: pl.DataFrame, probe_annotation: pl.DataFrame
) -> None:
    mapper = VariantTissueMapper()
    links = mapper.map_methylation_markers(annotated_variants, probe_annotation)
    genes = set(links["gene_symbol"].to_list())
    assert "FOXO3" in genes
    assert "APOE" in genes
    assert "CETP" not in genes
    assert "cg0001" in links["IlmnID"].to_list()


def test_map_variants_to_tissues_full(
    annotated_variants: pl.DataFrame,
    eqtl_table: pl.DataFrame,
    probe_annotation: pl.DataFrame,
) -> None:
    mapper = VariantTissueMapper(target_tissues=DEFAULT_TARGET_TISSUES)
    result = mapper.map_variants_to_tissues(
        annotated_variants,
        eqtl_table,
        alphagenome=pl.DataFrame(
            {"rsid": ["rs1"], "alphagenome_diff": [0.2]}
        ),
        probe_annotation=probe_annotation,
    )
    assert set(result) >= {"annotated", "eqtl_summary", "methylation_links"}
    assert result["annotated"].height == 3
    assert result["methylation_links"].height >= 2


def test_phenotype_channel_scores(annotated_variants: pl.DataFrame) -> None:
    integrator = PhenotypeIntegrator()
    vep = integrator.score_vep_impact(annotated_variants)
    assert vep.to_list() == [
        VEP_IMPACT_SCORES["HIGH"],
        VEP_IMPACT_SCORES["MODERATE"],
        VEP_IMPACT_SCORES["LOW"],
    ]
    ag = integrator.score_alphagenome(annotated_variants)
    assert ag[0] == pytest.approx(0.5)  # 25/50
    am = integrator.score_alphamissense(annotated_variants)
    assert am[0] == pytest.approx(0.8)
    assert am[2] == pytest.approx(0.0)


def test_compute_composite_risk_bounds(annotated_variants: pl.DataFrame) -> None:
    integrator = PhenotypeIntegrator(weights=DEFAULT_WEIGHTS)
    # Attach fake gtex summary columns.
    frame = annotated_variants.with_columns(
        pl.Series("gtex_n_eqtls", [10, 0, 2]),
        pl.Series("gtex_best_p_value", [1e-10, None, 1e-3]),
    )
    scored = integrator.compute_composite_risk(frame)
    assert "composite_risk" in scored.columns
    risks = scored["composite_risk"].to_list()
    assert all(0.0 <= r <= 1.0 for r in risks)
    assert risks[0] > risks[2]


def test_integrate_sample_profiles() -> None:
    integrator = PhenotypeIntegrator()
    variant_risks = pl.DataFrame(
        {
            "rsid": ["rs1", "rs2"],
            "composite_risk": [0.8, 0.2],
        }
    )
    samples = pl.DataFrame(
        {
            "sample_id": ["S1", "S1", "S2"],
            "rsid": ["rs1", "rs2", "rs1"],
            "alt_dosage": [2, 1, 0],
            "cohort": ["A", "A", "B"],
        }
    )
    profiles = integrator.integrate_sample_profiles(variant_risks, samples)
    assert profiles.height == 2
    s1 = profiles.filter(pl.col("sample_id") == "S1").row(0, named=True)
    # dosage-weighted: (2*0.8 + 1*0.2) / (2+1) = 1.8/3 = 0.6
    assert s1["sample_risk"] == pytest.approx(0.6)
    assert s1["cohort"] == "A"
    s2 = profiles.filter(pl.col("sample_id") == "S2").row(0, named=True)
    assert s2["sample_risk"] == pytest.approx(0.0)


def test_integrate_sample_profiles_requires_composite_risk() -> None:
    integrator = PhenotypeIntegrator()
    with pytest.raises(ValueError, match="composite_risk"):
        integrator.integrate_sample_profiles(
            pl.DataFrame({"rsid": ["rs1"]}),
            pl.DataFrame(
                {"sample_id": ["S1"], "rsid": ["rs1"], "alt_dosage": [1]}
            ),
        )


def test_run_integrative_pipeline(
    annotated_variants: pl.DataFrame,
    eqtl_table: pl.DataFrame,
    probe_annotation: pl.DataFrame,
) -> None:
    samples = pl.DataFrame(
        {
            "sample_id": ["S1", "S1", "S2"],
            "rsid": ["rs1", "rs2", "rs3"],
            "alt_dosage": [1, 2, 1],
        }
    )
    result = run_integrative_pipeline(
        annotated_variants,
        eqtl_table,
        alphagenome=pl.DataFrame(
            {"rsid": ["rs1", "rs2"], "perc_change": [40.0, -10.0]}
        ),
        probe_annotation=probe_annotation,
        sample_phenotypes=samples,
    )
    assert "variant_risks" in result
    assert "sample_profiles" in result
    assert "methylation_links" in result
    assert "composite_risk" in result["variant_risks"].columns
    assert result["sample_profiles"].height == 2
    assert isinstance(result["mapper"], VariantTissueMapper)
    assert isinstance(result["integrator"], PhenotypeIntegrator)


def test_package_exports() -> None:
    import rogen_aging.integrative as integrative

    assert "run_integrative_pipeline" in integrative.__all__
    assert callable(integrative.run_integrative_pipeline)


def test_summarize_eqtls_empty_after_filter() -> None:
    mapper = VariantTissueMapper(target_tissues=("Whole_Blood",))
    eqtls = pl.DataFrame(
        {
            "rsid": ["rs1"],
            "tissue": ["Lung"],
            "nes": [0.5],
            "p_value": [1e-5],
        }
    )
    summary = mapper.summarize_eqtl_profiles(eqtls)
    assert summary.is_empty()
    assert "gtex_n_eqtls" in summary.columns


def test_filter_eqtls_missing_tissue_raises() -> None:
    mapper = VariantTissueMapper()
    with pytest.raises(ValueError, match="missing required column: tissue"):
        mapper.filter_eqtls_to_target_tissues(pl.DataFrame({"rsid": ["rs1"]}))


def test_summarize_eqtls_missing_cols_raises() -> None:
    mapper = VariantTissueMapper()
    with pytest.raises(ValueError, match="missing required columns"):
        mapper.summarize_eqtl_profiles(pl.DataFrame({"rsid": ["rs1"], "tissue": ["Whole_Blood"]}))


def test_join_alphagenome_requires_join_key() -> None:
    mapper = VariantTissueMapper()
    variants = pl.DataFrame(
        {"chrom": ["1"], "pos": [1], "ref": ["A"], "alt": ["G"]}
    )
    scores = pl.DataFrame({"alphagenome_diff": [0.1]})
    with pytest.raises(ValueError, match="variant_key or rsid"):
        mapper.join_alphagenome_scores(variants, scores)


def test_map_methylation_missing_gene_raises(
    annotated_variants: pl.DataFrame, probe_annotation: pl.DataFrame
) -> None:
    mapper = VariantTissueMapper()
    bare = annotated_variants.select(["chrom", "pos", "ref", "alt"])
    with pytest.raises(ValueError, match="missing gene column"):
        mapper.map_methylation_markers(bare, probe_annotation)


def test_score_channels_without_optional_columns() -> None:
    integrator = PhenotypeIntegrator()
    bare = pl.DataFrame({"chrom": ["1"], "pos": [1], "ref": ["A"], "alt": ["G"]})
    assert integrator.score_vep_impact(bare).to_list() == [0.0]
    assert integrator.score_alphagenome(bare).to_list() == [0.0]
    assert integrator.score_alphamissense(bare).to_list() == [0.0]
    assert integrator.score_gtex_eqtl(bare).to_list() == [0.0]
    assert integrator.score_epigenetic(bare).to_list() == [0.0]


def test_score_alphagenome_perc_change_and_diff() -> None:
    integrator = PhenotypeIntegrator()
    perc = pl.DataFrame({"alphagenome_perc_change": [25.0, -100.0]})
    scored = integrator.score_alphagenome(perc).to_list()
    assert scored[0] == pytest.approx(0.5)
    assert scored[1] == pytest.approx(1.0)

    diff = pl.DataFrame({"alphagenome_diff": [0.4, -2.0]})
    scored_diff = integrator.score_alphagenome(diff).to_list()
    assert scored_diff[0] == pytest.approx(0.4)
    assert scored_diff[1] == pytest.approx(1.0)


def test_score_gtex_hit_count_only() -> None:
    integrator = PhenotypeIntegrator()
    frame = pl.DataFrame({"gtex_n_eqtls": [0, 5, 100]})
    scores = integrator.score_gtex_eqtl(frame).to_list()
    assert scores[0] == pytest.approx(0.0)
    assert scores[1] == pytest.approx(0.5)
    assert 0.0 < scores[2] <= 1.0


def test_compute_composite_risk_zero_weights_raises(
    annotated_variants: pl.DataFrame,
) -> None:
    integrator = PhenotypeIntegrator(
        weights={
            "vep_impact": 0.0,
            "alphagenome": 0.0,
            "alphamissense": 0.0,
            "gtex_eqtl": 0.0,
            "epigenetic": 0.0,
        }
    )
    with pytest.raises(ValueError, match="positive value"):
        integrator.compute_composite_risk(annotated_variants)


def test_build_risk_profile_without_samples(annotated_variants: pl.DataFrame) -> None:
    integrator = PhenotypeIntegrator()
    result = integrator.build_risk_profile(annotated_variants)
    assert set(result) == {"variant_risks"}
    assert "composite_risk" in result["variant_risks"].columns


def test_integrate_sample_profiles_missing_dosage_raises() -> None:
    integrator = PhenotypeIntegrator()
    with pytest.raises(ValueError, match="alt_dosage"):
        integrator.integrate_sample_profiles(
            pl.DataFrame({"rsid": ["rs1"], "composite_risk": [0.5]}),
            pl.DataFrame({"sample_id": ["S1"], "rsid": ["rs1"]}),
        )


def test_run_integrative_pipeline_without_optional_tables(
    annotated_variants: pl.DataFrame, eqtl_table: pl.DataFrame
) -> None:
    result = run_integrative_pipeline(annotated_variants, eqtl_table)
    assert "variant_risks" in result
    assert "sample_profiles" not in result
    assert "methylation_links" not in result
    assert result["variant_risks"].height == annotated_variants.height


def test_normalize_rsid_edge_tokens() -> None:
    assert normalize_rsid("  ") is None
    assert normalize_rsid("nan") is None
    assert normalize_rsid("None") is None
    assert normalize_rsid("rs0") == "rs0"


def test_variant_key_multiallelic_uses_first_alt() -> None:
    assert variant_key("chrX", 10, "a", " t , g ") == "X:10:A:T"
