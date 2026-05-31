"""Tests for ``rogen_aging.ukb.mock_rap`` UKB-RAP synthetic layout."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from rogen_aging.ukb.mock_rap import (
    PHENOTYPE_V2_FIELDS,
    generate_ukb_rap_mock,
    load_snp_manifest,
)

LA_SNP_COUNT = 70


@pytest.fixture
def la_snp_manifest_70(tmp_path: Path) -> Path:
    """Minimal 70-row manifest with GRCh38 coordinates on chr1–chr22."""
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
def mock_ukb_rap_dir(tmp_path: Path, la_snp_manifest_70: Path) -> Path:
    generate_ukb_rap_mock(
        n_samples=16,
        snp_manifest=la_snp_manifest_70,
        output_dir=tmp_path / "mock_ukb_rap",
        seed=11,
    )
    return tmp_path / "mock_ukb_rap"


def test_manifest_loads_seventy_snps(la_snp_manifest_70: Path) -> None:
    manifest = load_snp_manifest(la_snp_manifest_70)
    assert len(manifest) == LA_SNP_COUNT


def test_phenotype_and_vcf_share_eid_set(mock_ukb_rap_dir: Path) -> None:
    phenotype_path = mock_ukb_rap_dir / "phenotypes" / "ukb_phenotypes.csv"
    vcf_path = mock_ukb_rap_dir / "genotypes" / "ukb_la_snps.vcf"

    pheno = pl.read_csv(phenotype_path, comment_prefix="#")
    vcf_lines = vcf_path.read_text(encoding="ascii").splitlines()
    header_line = next(line for line in vcf_lines if line.startswith("#CHROM"))
    vcf_sample_ids = header_line.split("\t")[9:]

    pheno_eids = set(pheno["eid"].to_list())
    vcf_eids = set(vcf_sample_ids)
    assert pheno_eids == vcf_eids
    assert len(pheno_eids) == 16


def test_vcf_has_exactly_seventy_variants(mock_ukb_rap_dir: Path) -> None:
    vcf_path = mock_ukb_rap_dir / "genotypes" / "ukb_la_snps.vcf"
    lines = vcf_path.read_text(encoding="ascii").splitlines()
    data_lines = [line for line in lines if line and not line.startswith("#")]
    assert len(data_lines) == LA_SNP_COUNT


def test_phenotype_v2_columns_present(mock_ukb_rap_dir: Path) -> None:
    phenotype_path = mock_ukb_rap_dir / "phenotypes" / "ukb_phenotypes.csv"
    pheno = pl.read_csv(phenotype_path, comment_prefix="#")
    assert "eid" in pheno.columns
    for field in PHENOTYPE_V2_FIELDS:
        assert field in pheno.columns, f"missing v2 phenotype field: {field}"
