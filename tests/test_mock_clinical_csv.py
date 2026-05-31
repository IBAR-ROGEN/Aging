"""Pytest checks for synthetic cohort tabular output from the mock clinical CSV generator."""

from __future__ import annotations

import polars as pl
import pytest

from rogen_aging.ukb.mock_clinical import DUMMY_SNP_IDS, generate_synthetic_ukb_data

REQUIRED_CLINICAL_COLUMNS: tuple[str, ...] = (
    "Sample_ID",
    "Age",
    "Sex",
    "BMI",
    "AD_diagnosis",
    "EAA",
)


@pytest.fixture
def synthetic_clinical_frame() -> pl.DataFrame:
    """Small synthetic cohort in a biologically wide adult age range."""
    df = generate_synthetic_ukb_data(
        n_samples=32,
        min_age=18,
        max_age=110,
        seed=7,
    )
    return pl.from_pandas(df)


def test_mock_clinical_required_columns(synthetic_clinical_frame: pl.DataFrame) -> None:
    cols = set(synthetic_clinical_frame.columns)
    for name in REQUIRED_CLINICAL_COLUMNS:
        assert name in cols, f"missing required column: {name}"
    for snp in DUMMY_SNP_IDS:
        assert snp in cols, f"missing SNP column: {snp}"


def test_mock_clinical_plausible_ranges(synthetic_clinical_frame: pl.DataFrame) -> None:
    age = synthetic_clinical_frame["Age"]
    assert age.min() >= 18
    assert age.max() <= 110

    bmi = synthetic_clinical_frame["BMI"]
    assert bmi.min() >= 15.0
    assert bmi.max() <= 50.0

    sex = synthetic_clinical_frame["Sex"]
    assert set(sex.unique().to_list()) <= {0, 1}

    ad = synthetic_clinical_frame["AD_diagnosis"]
    assert set(ad.unique().to_list()) <= {0, 1}

    for snp in DUMMY_SNP_IDS:
        g = synthetic_clinical_frame[snp]
        assert g.min() >= 0
        assert g.max() <= 2
        assert g.dtype.is_integer()

    assert synthetic_clinical_frame["Sample_ID"].str.starts_with("MOCK_").all()
    assert synthetic_clinical_frame.height == 32
