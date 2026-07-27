"""Tests for Activity 2.1.11.1 synthetic UKB integrative join + association."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from rogen_aging.ukb.mock_rap import generate_ukb_rap_mock
from rogen_aging.ukb_integration.ukb_joiner import (
    AD_OUT,
    LA_SNP_ASSOC_COLUMNS,
    PARENTAL_LONGEVITY_OUT,
    join_phenotypes_genotypes,
    load_genotype_matrix_from_vcf,
    load_phenotype_table,
    run_association_scan,
    run_integration_pipeline,
)

LA_SNP_COUNT = 70


@pytest.fixture
def la_snp_manifest_70(tmp_path: Path) -> Path:
    rows: list[dict[str, object]] = []
    for idx in range(LA_SNP_COUNT):
        chrom_num = (idx % 22) + 1
        rows.append(
            {
                "Gene": f"GENE_{idx + 1:03d}",
                "SNP_rsID": f"rs_mock_{idx + 1:04d}",
                "Chromosome": str(chrom_num),
                "Position_GRCh38": 10_000 + idx * 1_000,
            }
        )
    path = tmp_path / "mock_la_snp_manifest_70.csv"
    pl.DataFrame(rows).write_csv(path)
    return path


@pytest.fixture
def mock_ukb_rap_paths(tmp_path: Path, la_snp_manifest_70: Path) -> tuple[Path, Path, Path]:
    root = tmp_path / "mock_ukb_rap"
    generate_ukb_rap_mock(
        n_samples=24,
        snp_manifest=la_snp_manifest_70,
        output_dir=root,
        seed=7,
    )
    pheno = root / "phenotypes" / "ukb_phenotypes.csv"
    vcf = root / "genotypes" / "ukb_la_snps.vcf"
    return pheno, vcf, tmp_path / "out"


def test_joiner_one_row_per_eid(mock_ukb_rap_paths: tuple[Path, Path, Path]) -> None:
    pheno_path, vcf_path, _ = mock_ukb_rap_paths
    phenotypes = load_phenotype_table(pheno_path)
    genotypes = load_genotype_matrix_from_vcf(vcf_path)
    joined = join_phenotypes_genotypes(phenotypes, genotypes)

    assert joined.height == phenotypes.height
    assert joined["eid"].n_unique() == joined.height
    assert len([c for c in joined.columns if c.startswith("rs_mock_")]) == LA_SNP_COUNT


def test_association_outputs_have_seventy_rows(
    mock_ukb_rap_paths: tuple[Path, Path, Path],
) -> None:
    pheno_path, vcf_path, out_dir = mock_ukb_rap_paths
    joined, parental, ad = run_integration_pipeline(pheno_path, vcf_path, out_dir)

    assert parental.height == LA_SNP_COUNT
    assert ad.height == LA_SNP_COUNT
    for frame in (parental, ad):
        assert list(frame.columns) == list(LA_SNP_ASSOC_COLUMNS)

    parental_path = out_dir / PARENTAL_LONGEVITY_OUT
    ad_path = out_dir / AD_OUT
    assert parental_path.is_file()
    assert ad_path.is_file()

    reloaded_parental = pl.read_csv(parental_path, comment_prefix="#")
    assert reloaded_parental.height == LA_SNP_COUNT

    snp_cols = [c for c in joined.columns if c.startswith("rs_mock_")]
    scan = run_association_scan(joined, phenotype_col="parental_longevity")
    assert scan.height == len(snp_cols) == LA_SNP_COUNT


def test_alt_dosage_from_gt_type_cyvcf2_encoding() -> None:
    from rogen_aging.ukb_integration.ukb_joiner import _alt_dosage_from_gt_type

    # cyvcf2: 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
    assert _alt_dosage_from_gt_type(0) == 0
    assert _alt_dosage_from_gt_type(1) == 1
    assert _alt_dosage_from_gt_type(3) == 2
    assert _alt_dosage_from_gt_type(2) is None


def test_load_genotype_matrix_maps_hom_alt(tmp_path: Path) -> None:
    vcf_text = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
1	100	rsTEST	A	G	.	.	.	GT	0/0	0/1	1/1	./.
"""
    path = tmp_path / "gt.vcf"
    path.write_text(vcf_text, encoding="utf-8")
    frame = load_genotype_matrix_from_vcf(path)
    assert frame["rsTEST"].to_list() == [0, 1, 2, None]


def test_contingency_drops_missing_dosages() -> None:
    import numpy as np

    from rogen_aging.ukb_integration.ukb_joiner import genotype_phenotype_contingency

    g = np.array([0.0, 1.0, 2.0, np.nan])
    y = np.array([0, 1, 0, 1])
    table = genotype_phenotype_contingency(g, y)
    # Missing dosage must not become HOM_REF under outcome=1.
    assert table.tolist() == [[1, 0, 1], [0, 1, 0]]
    assert int(table.sum()) == 3