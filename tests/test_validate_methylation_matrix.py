"""Tests for methylation matrix preflight validation."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from rogen_aging.clock.validate_matrix import (
    MissingValuePolicy,
    validate_methylation_matrix,
    write_validation_outputs,
)


def _toy_tables(
    tmp_path: Path,
    *,
    with_missing: bool = False,
    with_oor: bool = False,
    drop_meta_sample: bool = False,
    drop_matrix_cpg: bool = False,
) -> tuple[Path, Path, Path]:
    sample_ids = ["S1", "S2", "S3"]
    cpgs = ["cg0001", "cg0002", "cg0003"]
    rows = [
        [0.21, 0.55, 0.77],
        [0.33, 0.44, 0.66],
        [0.12, 0.58, 0.91],
    ]
    if with_missing:
        rows[1][1] = float("nan")
    if with_oor:
        rows[0][0] = 1.2
    if drop_matrix_cpg:
        cpgs = cpgs[:-1]
        rows = [r[:-1] for r in rows]

    matrix = pd.DataFrame(rows, columns=cpgs)
    matrix.insert(0, "sample_id", sample_ids)
    meta_ids = sample_ids[:-1] if drop_meta_sample else sample_ids
    meta = pd.DataFrame(
        {
            "sample_id": meta_ids,
            "chronological_age": [40.0, 55.0, 62.0][: len(meta_ids)],
        }
    )
    expected = tmp_path / "expected_cpgs.txt"
    expected.write_text("\n".join(["cg0001", "cg0002", "cg0003"]) + "\n", encoding="utf-8")
    matrix_path = tmp_path / "matrix.csv"
    meta_path = tmp_path / "meta.csv"
    matrix.to_csv(matrix_path, index=False)
    meta.to_csv(meta_path, index=False)
    return matrix_path, meta_path, expected


def test_validate_clean_matrix_passes(tmp_path: Path) -> None:
    matrix_path, meta_path, expected = _toy_tables(tmp_path)
    report, cleaned = validate_methylation_matrix(
        matrix_path,
        meta_path,
        expected_cpgs=expected,
        missing_policy=MissingValuePolicy.FAIL,
    )
    assert report.passed
    assert report.n_samples == 3
    assert report.n_cpgs == 3
    assert report.n_missing_beta_values == 0
    assert report.n_out_of_range_beta_values == 0
    assert list(cleaned.columns) == ["sample_id", "cg0001", "cg0002", "cg0003"]
    log = report.format_log()
    assert "PASSED" in log
    assert "Samples: 3" in log


def test_missing_values_impute_and_report(tmp_path: Path) -> None:
    matrix_path, meta_path, expected = _toy_tables(tmp_path, with_missing=True)
    report, cleaned = validate_methylation_matrix(
        matrix_path,
        meta_path,
        expected_cpgs=expected,
        missing_policy=MissingValuePolicy.IMPUTE_COLUMN_MEAN,
    )
    assert report.passed  # missingness is warning under impute
    assert report.n_missing_beta_values == 1
    assert cleaned["cg0002"].isna().sum() == 0
    assert any(i.check == "missing_values" for i in report.issues)


def test_missing_values_fail_policy(tmp_path: Path) -> None:
    matrix_path, meta_path, expected = _toy_tables(tmp_path, with_missing=True)
    report, _ = validate_methylation_matrix(
        matrix_path,
        meta_path,
        expected_cpgs=expected,
        missing_policy=MissingValuePolicy.FAIL,
    )
    assert not report.passed
    assert any(i.check == "missing_values" and i.severity == "error" for i in report.issues)


def test_beta_out_of_range_fails(tmp_path: Path) -> None:
    matrix_path, meta_path, expected = _toy_tables(tmp_path, with_oor=True)
    report, _ = validate_methylation_matrix(
        matrix_path,
        meta_path,
        expected_cpgs=expected,
    )
    assert not report.passed
    assert report.n_out_of_range_beta_values == 1
    assert any(i.check == "beta_range" for i in report.issues)


def test_sample_id_manifest_mismatch(tmp_path: Path) -> None:
    matrix_path, meta_path, expected = _toy_tables(tmp_path, drop_meta_sample=True)
    report, _ = validate_methylation_matrix(
        matrix_path,
        meta_path,
        expected_cpgs=expected,
    )
    assert not report.passed
    assert report.n_samples_in_matrix_only == 1
    assert any(i.check == "sample_ids" for i in report.issues)


def test_missing_expected_cpg(tmp_path: Path) -> None:
    matrix_path, meta_path, expected = _toy_tables(tmp_path, drop_matrix_cpg=True)
    report, _ = validate_methylation_matrix(
        matrix_path,
        meta_path,
        expected_cpgs=expected,
    )
    assert not report.passed
    assert report.n_missing_expected_cpgs == 1
    assert any(i.check == "expected_cpgs" for i in report.issues)


def test_write_outputs(tmp_path: Path) -> None:
    matrix_path, meta_path, expected = _toy_tables(tmp_path)
    report, cleaned = validate_methylation_matrix(
        matrix_path,
        meta_path,
        expected_cpgs=expected,
    )
    log_path = tmp_path / "preflight.log"
    json_path = tmp_path / "preflight.json"
    cleaned_path = tmp_path / "cleaned.csv"
    write_validation_outputs(
        report,
        cleaned,
        log_path=log_path,
        report_json_path=json_path,
        cleaned_matrix_path=cleaned_path,
    )
    assert log_path.is_file()
    assert "PASSED" in log_path.read_text(encoding="utf-8")
    assert json_path.is_file()
    assert cleaned_path.is_file()
    assert "sample_id" in pd.read_csv(cleaned_path).columns


def test_boundary_beta_fails_strict(tmp_path: Path) -> None:
    matrix_path, meta_path, expected = _toy_tables(tmp_path)
    matrix = pd.read_csv(matrix_path)
    matrix.loc[0, "cg0001"] = 0.0
    matrix.to_csv(matrix_path, index=False)
    report, _ = validate_methylation_matrix(
        matrix_path,
        meta_path,
        expected_cpgs=expected,
        require_strict_beta=True,
    )
    assert not report.passed
    assert report.n_out_of_range_beta_values == 1
