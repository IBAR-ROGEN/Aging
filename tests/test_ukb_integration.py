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
    JoinDropRateError,
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


def test_joiner_one_row_per_eid(
    mock_ukb_rap_paths: tuple[Path, Path, Path],
    tmp_path: Path,
) -> None:
    pheno_path, vcf_path, _ = mock_ukb_rap_paths
    phenotypes = load_phenotype_table(pheno_path)
    genotypes = load_genotype_matrix_from_vcf(vcf_path)
    audit_log = tmp_path / "ukb_integration_audit.log"
    joined = join_phenotypes_genotypes(phenotypes, genotypes, audit_log=audit_log)

    assert joined.height == phenotypes.height
    assert joined["eid"].n_unique() == joined.height
    assert len([c for c in joined.columns if c.startswith("rs_mock_")]) == LA_SNP_COUNT
    assert audit_log.is_file()
    audit_text = audit_log.read_text(encoding="utf-8")
    assert "drop_rate=0.000000" in audit_text
    assert "halted=False" in audit_text


def test_association_outputs_have_seventy_rows(
    mock_ukb_rap_paths: tuple[Path, Path, Path],
    tmp_path: Path,
) -> None:
    pheno_path, vcf_path, out_dir = mock_ukb_rap_paths
    audit_log = tmp_path / "ukb_integration_audit.log"
    joined, parental, ad = run_integration_pipeline(
        pheno_path, vcf_path, out_dir, audit_log=audit_log
    )

    assert parental.height == LA_SNP_COUNT
    assert ad.height == LA_SNP_COUNT
    for frame in (parental, ad):
        assert list(frame.columns) == list(LA_SNP_ASSOC_COLUMNS)

    parental_path = out_dir / PARENTAL_LONGEVITY_OUT
    ad_path = out_dir / AD_OUT
    assert parental_path.is_file()
    assert ad_path.is_file()
    assert audit_log.is_file()

    reloaded_parental = pl.read_csv(parental_path, comment_prefix="#")
    assert reloaded_parental.height == LA_SNP_COUNT

    snp_cols = [c for c in joined.columns if c.startswith("rs_mock_")]
    scan = run_association_scan(joined, phenotype_col="parental_longevity")
    assert scan.height == len(snp_cols) == LA_SNP_COUNT


def test_eid_schema_mismatch_normalized_and_audited(tmp_path: Path) -> None:
    phenotypes = pl.DataFrame(
        {
            "eid": [1001, 1002, 1003],
            "parental_longevity": [0, 1, 0],
            "ad_diagnosis_code": ["", "G30.1", ""],
        }
    )
    genotypes = pl.DataFrame(
        {
            "eid": ["1001", "1002", "1003"],
            "rs_mock_0001": [0, 1, 2],
        }
    )
    audit_log = tmp_path / "ukb_integration_audit.log"
    joined = join_phenotypes_genotypes(phenotypes, genotypes, audit_log=audit_log)

    assert joined.height == 3
    assert joined["eid"].to_list() == ["1001", "1002", "1003"]
    audit_text = audit_log.read_text(encoding="utf-8")
    assert "eid_schema_mismatch=True" in audit_text
    assert "eid_schema_normalized_to_utf8=True" in audit_text
    assert "drop_rate=0.000000" in audit_text


def test_join_halts_when_drop_rate_exceeds_one_percent(tmp_path: Path) -> None:
    phenotypes = pl.DataFrame(
        {
            "eid": [f"E{i}" for i in range(100)],
            "parental_longevity": [0] * 100,
            "ad_diagnosis_code": [""] * 100,
        }
    )
    # Only 97 shared IDs → 3 pheno-only + 0 geno-only over union 100 → 3% drop.
    genotypes = pl.DataFrame(
        {
            "eid": [f"E{i}" for i in range(97)],
            "rs_mock_0001": [0] * 97,
        }
    )
    audit_log = tmp_path / "ukb_integration_audit.log"
    with pytest.raises(JoinDropRateError, match="drop rate"):
        join_phenotypes_genotypes(phenotypes, genotypes, audit_log=audit_log)

    audit_text = audit_log.read_text(encoding="utf-8")
    assert "halted=True" in audit_text
    assert "pheno_only_dropped=3" in audit_text


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