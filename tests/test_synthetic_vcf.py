"""Tests for ``scripts/generate_synthetic_romanian_vcf.py`` synthetic VCF output."""

from __future__ import annotations

from pathlib import Path

import pytest

from generate_synthetic_romanian_vcf import main as write_synthetic_vcf


@pytest.fixture
def tiny_vcf_path(tmp_path: Path) -> Path:
    out = tmp_path / "synthetic_minimal.vcf"
    write_synthetic_vcf(
        [
            "--samples",
            "4",
            "--variants",
            "3",
            "--output",
            str(out),
            "--seed",
            "99",
        ]
    )
    return out


def test_synthetic_vcf_headers_and_format(tiny_vcf_path: Path) -> None:
    text = tiny_vcf_path.read_text(encoding="ascii")
    assert "##fileformat=VCFv4.2\n" in text
    assert "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" in text


def test_synthetic_vcf_column_count_matches_samples(tiny_vcf_path: Path) -> None:
    lines = tiny_vcf_path.read_text(encoding="ascii").splitlines()
    header_line = next(line for line in lines if line.startswith("#CHROM"))
    assert header_line.startswith("#CHROM\t")

    fixed_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    header_parts = header_line.split("\t")
    assert header_parts[: len(fixed_cols)] == fixed_cols
    sample_cols = header_parts[len(fixed_cols) :]
    n_samples = len(sample_cols)
    assert n_samples == 4

    data_lines = [ln for ln in lines if ln and not ln.startswith("#")]
    assert len(data_lines) == 3
    for row in data_lines:
        parts = row.split("\t")
        assert len(parts) == len(fixed_cols) + n_samples
