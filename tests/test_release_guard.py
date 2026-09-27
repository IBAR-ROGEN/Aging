"""Release-guard checks on small tables built in the test."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from rogen_aging.cli.release_check import app
from rogen_aging.config import default_config_path, load_config
from rogen_aging.privacy.release_guard import run_release_guard

_RUNNER = CliRunner()


def _check(frame: pl.DataFrame, kind: str, name: str, *, n_samples: int | None = None):
    result = run_release_guard(
        frame,
        table_kind=kind,
        min_cell_count=5,
        n_samples=n_samples,
    )
    return next(item for item in result.checks if item.check == name), result


def test_config_min_cell_count_is_project_default() -> None:
    text = default_config_path().read_text(encoding="utf-8")
    assert "project default, to be confirmed with the data custodian" in text
    cfg = load_config()
    assert int(cfg.release_guard.min_cell_count) == 5


def test_individual_level_column_fails() -> None:
    frame = pl.DataFrame({"sample_id": ["a", "b"], "gene": ["APOE", "CETP"]})
    check, result = _check(frame, "strata", "individual_level")
    assert check.status == "FAIL"
    assert "sample_id" in check.reason
    assert result.release_blocked
    assert result.sanitized is None


def test_pseudonym_in_free_text_fails() -> None:
    frame = pl.DataFrame(
        {
            "gene": ["APOE"],
            "note": ["mentioned RO-abcdef012345 in the free-text comment"],
        }
    )
    check, result = _check(frame, "strata", "individual_level")
    assert check.status == "FAIL"
    assert "note" in check.reason
    assert "RO-" in check.reason
    assert result.release_blocked
    assert result.sanitized is None


def test_one_row_per_sample_fails() -> None:
    frame = pl.DataFrame({"gene": ["APOE", "CETP", "TLR4"], "n_carriers": [10, 12, 8]})
    check, result = _check(frame, "strata", "individual_level", n_samples=3)
    assert check.status == "FAIL"
    assert "one row per sample" in check.reason
    assert result.release_blocked


def test_direct_identifier_columns_fail() -> None:
    frame = pl.DataFrame(
        {
            "date_of_birth": ["1990-01-15"],
            "postcode": ["SW1A 1AA"],
            "full_name": ["Ada Lovelace"],
            "gene": ["APOE"],
        }
    )
    check, result = _check(frame, "strata", "direct_identifiers")
    assert check.status == "FAIL"
    assert "date of birth" in check.reason
    assert "postcode" in check.reason
    assert "full_name" in check.reason
    assert result.sanitized is None


def test_allele_complementary_suppression() -> None:
    frame = pl.DataFrame(
        {
            "variant": ["small_ac", "small_complement", "kept"],
            "AC": [2, 40, 10],
            "AN": [100, 43, 100],
            "AF": [0.02, 0.930232, 0.123456],
        }
    )
    check, result = _check(frame, "allele_counts", "allele_count_cells")
    assert check.status == "PASS"
    assert "suppressed 2 row(s)" in check.reason
    assert "AN - AC < 5" in check.reason
    assert result.sanitized is not None
    assert result.sanitized["variant"].to_list() == ["kept"]
    assert result.sanitized["AC"].to_list() == [10]
    assert result.sanitized["AF"][0] == pytest.approx(0.123)


def test_summary_kind_leaves_frequencies_unchanged() -> None:
    frame = pl.DataFrame({"rsid": ["rs1"], "maf_nfe": [0.01], "power": [0.2]})
    check, result = _check(frame, "summary", "strata_cells")
    assert check.status == "PASS"
    assert "not applied: table kind is summary" in check.reason
    assert result.sanitized is not None
    assert result.sanitized["maf_nfe"].to_list() == [0.01]
    assert result.sanitized["power"].to_list() == [0.2]
    assert result.release_blocked is False


def test_strata_complementary_suppression() -> None:
    frame = pl.DataFrame(
        {
            "phenotype": ["need_extra", "already_two"],
            "group_a": [4, 1],
            "group_b": [10, 2],
            "group_c": [30, 40],
        }
    )
    check, result = _check(frame, "strata", "strata_cells")
    assert check.status == "PASS"
    assert "complementary suppression hid 1" in check.reason
    assert result.sanitized is not None
    sanitized = result.sanitized
    assert sanitized["group_a"].null_count() == 2
    assert sanitized["group_b"].null_count() == 2
    assert sanitized["group_c"].to_list() == [30, 40]
    assert sanitized["phenotype"].to_list() == ["need_extra", "already_two"]


def test_write_sanitized_drops_complementary_small_counts(tmp_path: Path) -> None:
    table = tmp_path / "ac.csv"
    report = tmp_path / "report.json"
    sanitized = tmp_path / "clean.csv"
    pl.DataFrame(
        {
            "variant": ["small_complement", "kept"],
            "AC": [40, 10],
            "AN": [43, 100],
            "AF": [0.930232, 0.123456],
        }
    ).write_csv(table)
    result = _RUNNER.invoke(
        app,
        [
            "--table",
            str(table),
            "--table-kind",
            "allele_counts",
            "--report",
            str(report),
            "--write-sanitized",
            str(sanitized),
        ],
    )
    assert result.exit_code == 0, result.output
    assert sanitized.is_file()
    kept = pl.read_csv(sanitized)
    assert kept["variant"].to_list() == ["kept"]
    assert kept["AF"][0] == pytest.approx(0.123)


def test_identifier_failure_does_not_write_sanitized(tmp_path: Path) -> None:
    table = tmp_path / "notes.csv"
    report = tmp_path / "report.json"
    sanitized = tmp_path / "clean.csv"
    pl.DataFrame({"note": ["RO-abcdef012345 slipped into text"], "n": [20]}).write_csv(table)
    result = _RUNNER.invoke(
        app,
        [
            "--table",
            str(table),
            "--table-kind",
            "strata",
            "--report",
            str(report),
            "--write-sanitized",
            str(sanitized),
        ],
    )
    assert result.exit_code == 1
    assert report.is_file()
    assert not sanitized.exists()
    payload = json.loads(report.read_text(encoding="utf-8"))
    assert payload["release_blocked"] is True
    assert payload["sanitized_written"] is False
    assert [item["check"] for item in payload["checks"]] == [
        "individual_level",
        "direct_identifiers",
        "allele_count_cells",
        "strata_cells",
    ]
