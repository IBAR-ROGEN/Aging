"""Unit tests for GEO series-matrix sex / ethnicity parsing."""

from __future__ import annotations

from pathlib import Path

from rogen_aging.clock.external_data import (
    load_series_phenotypes,
    normalize_sex_label,
    parse_sample_phenotypes,
)


def _tiny_series_matrix(path: Path) -> None:
    lines = [
        "!Sample_geo_accession\tGSM1\tGSM2\tGSM3",
        '!Sample_characteristics_ch1\t"gender: Male"\t"gender: F"\t"sex: female"',
        '!Sample_characteristics_ch1\t"age: 72"\t"age (y): 55"\t"age: 40.5"',
        '!Sample_characteristics_ch1\t"ethnicity: Caucasian - European"\t'
        '"ethnicity: African American"\t"tissue: whole blood"',
        "!series_matrix_table_begin",
        "ID_REF\tGSM1\tGSM2\tGSM3",
        "cg00000029\t0.1\t0.2\t0.3",
        "!series_matrix_table_end",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_normalize_sex_label() -> None:
    assert normalize_sex_label("Male") == "M"
    assert normalize_sex_label("f") == "F"
    assert normalize_sex_label("unknown") is None


def test_parse_sample_phenotypes_from_fixture(tmp_path: Path) -> None:
    matrix = tmp_path / "tiny_series_matrix.txt"
    _tiny_series_matrix(matrix)
    frame = load_series_phenotypes(matrix)
    assert list(frame["sample_id"]) == ["GSM1", "GSM2", "GSM3"]
    assert frame["chronological_age"].tolist() == [72.0, 55.0, 40.5]
    assert frame["sex"].tolist() == ["M", "F", "F"]
    assert frame["female"].tolist() == [0, 1, 1]
    assert frame["ethnicity"].tolist() == [
        "Caucasian - European",
        "African American",
        None,
    ]


def test_parse_sample_phenotypes_direct() -> None:
    sample_rows = {
        "!Sample_geo_accession": ["A", "B"],
        "!Sample_characteristics_ch1": [
            ["gender: Female", "gender: Male"],
            ["age: 10", "age: 20"],
        ],
    }
    frame = parse_sample_phenotypes(sample_rows)
    assert frame["sex"].tolist() == ["F", "M"]
    assert frame["female"].tolist() == [1, 0]
    assert frame["ethnicity"].tolist() == [None, None]
