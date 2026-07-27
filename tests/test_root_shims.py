"""Smoke tests for root deprecation shims and canonical script paths."""

from __future__ import annotations

from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent

_SHIMS = (
    ("plot_af_comparison.py", "scripts/figures/plot_af_comparison.py"),
    ("plot_clock_eval.py", "scripts/figures/plot_clock_eval.py"),
    ("plot_clock_validation.py", "scripts/figures/plot_clock_validation.py"),
    ("plot_consequence_summary.py", "scripts/figures/plot_consequence_summary.py"),
    ("evaluate_methylation_clock.py", "scripts/clock/evaluate_methylation_clock.py"),
    ("reconcile_and_generate_figures.py", "scripts/figures/reconcile_and_generate_figures.py"),
    ("run_july_annotation_pipeline.py", "scripts/ukb/run_july_annotation_pipeline.py"),
    ("annotate_la_snps_vep.py", "scripts/ukb/annotate_la_snps_vep.py"),
    ("annotate_la_snps_gtex.py", "scripts/ukb/annotate_la_snps_gtex.py"),
)


@pytest.mark.parametrize(("shim", "canonical"), _SHIMS)
def test_root_shim_forwards_to_canonical(shim: str, canonical: str) -> None:
    shim_path = _REPO_ROOT / shim
    canonical_path = _REPO_ROOT / canonical
    assert shim_path.is_file(), f"missing shim {shim_path}"
    assert canonical_path.is_file(), f"missing canonical {canonical_path}"
    text = shim_path.read_text(encoding="utf-8")
    assert "DeprecationWarning" in text
    assert "runpy.run_path" in text
    assert canonical.replace("\\", "/") in text.replace("\\", "/")


def test_canonical_figure_scripts_are_implementations() -> None:
    for relative in (
        "scripts/figures/plot_af_comparison.py",
        "scripts/figures/plot_clock_validation.py",
        "scripts/figures/plot_consequence_summary.py",
        "scripts/clock/evaluate_methylation_clock.py",
    ):
        path = _REPO_ROOT / relative
        text = path.read_text(encoding="utf-8")
        assert "runpy.run_path" not in text
        assert path.stat().st_size > 500
