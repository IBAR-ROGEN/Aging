"""Tests for 1KG vs gnomAD comparison summarization."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from rogen_aging.ukb.gnomad import format_summary_markdown, summarize_comparison


def test_summarize_comparison_counts_and_top_loci() -> None:
    comparison = pd.DataFrame(
        [
            {"rsID": "rs1", "AF_1kg": 0.10, "AF_gnomad_nfe": 0.11, "abs_diff": 0.01},
            {"rsID": "rs2", "AF_1kg": 0.20, "AF_gnomad_nfe": 0.30, "abs_diff": 0.10},
            {"rsID": "rs3", "AF_1kg": 0.05, "AF_gnomad_nfe": 0.10, "abs_diff": 0.05},
            {"rsID": "rs4", "AF_1kg": None, "AF_gnomad_nfe": 0.02, "abs_diff": None},
            {"rsID": "rs5", "AF_1kg": 0.03, "AF_gnomad_nfe": None, "abs_diff": None},
        ]
    )

    summary = summarize_comparison(comparison, diff_threshold=0.05)

    assert summary.total_snps == 5
    assert summary.resolved_both == 3
    assert summary.missing_1kg == 1
    assert summary.missing_gnomad == 1
    assert summary.concordant == 1
    assert summary.discordant == 2
    assert summary.mean_abs_diff_concordant == pytest.approx(0.01)
    assert summary.median_abs_diff_concordant == pytest.approx(0.01)
    assert summary.top_discordant.iloc[0]["rsID"] == "rs2"
    assert summary.top_discordant.iloc[1]["rsID"] == "rs3"


def test_format_summary_markdown_writes_table(tmp_path: Path) -> None:
    comparison = pd.DataFrame(
        [
            {"rsID": "rs1", "AF_1kg": 0.10, "AF_gnomad_nfe": 0.11},
            {"rsID": "rs2", "AF_1kg": 0.20, "AF_gnomad_nfe": 0.30},
        ]
    )
    summary = summarize_comparison(comparison)
    markdown = format_summary_markdown(summary)

    assert "Of 2 LA-SNPs in the comparison table" in markdown
    assert "| rsID | AF_1kg | AF_gnomad_nfe | abs_diff |" in markdown
    assert "rs2" in markdown

    output = tmp_path / "af_comparison_summary.md"
    output.write_text(markdown, encoding="utf-8")
    assert output.read_text(encoding="utf-8") == markdown
