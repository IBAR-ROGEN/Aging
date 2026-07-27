"""Smoke test for cluster ∩ LongevityMap overlap enrichment."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SCRIPT = _REPO_ROOT / "analysis" / "overlap_enrichment" / "run_overlap_enrichment.py"
_CLUSTER = _REPO_ROOT / "data" / "Supplementary Table 3.xlsx"
_LONGEVITY = _REPO_ROOT / "data" / "longevitymap.sqlite"
_SNPS = _REPO_ROOT / "results" / "snps_validated.csv"


@pytest.mark.skipif(not _SCRIPT.is_file(), reason="overlap enrichment script missing")
@pytest.mark.skipif(not _CLUSTER.is_file(), reason="cluster table missing")
@pytest.mark.skipif(not _LONGEVITY.is_file(), reason="longevity sqlite missing")
@pytest.mark.skipif(not _SNPS.is_file(), reason="snps_validated.csv missing")
def test_overlap_enrichment_run_match(tmp_path: Path) -> None:
    """End-to-end enrichment against local inputs must report MATCH (41 genes)."""
    out = tmp_path / "overlap_out"
    result = subprocess.run(
        [
            "uv",
            "run",
            "python",
            str(_SCRIPT),
            "run",
            "--offline",
            "--cluster-table",
            str(_CLUSTER),
            "--longevity-db",
            str(_LONGEVITY),
            "--snps-validated",
            str(_SNPS),
            "--output-dir",
            str(out),
            "--cache-dir",
            str(_REPO_ROOT / "results" / "cache"),
        ],
        cwd=_REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert "verdict=MATCH" in result.stdout
    report = (out / "overlap_enrichment.md").read_text(encoding="utf-8")
    assert "Verdict: **MATCH**" in report
    assert "Platform universe:" in report
    assert "SKIPPED" in report
    assert (out / "overlap_enrichment_stats.csv").is_file()
