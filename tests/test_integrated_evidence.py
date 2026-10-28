"""Unit tests for integrated evidence labels and power tiers."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import polars as pl

_REPO = Path(__file__).resolve().parent.parent
_SCRIPT = _REPO / "scripts" / "longevity" / "integrated_evidence.py"


def _load_script() -> ModuleType:
    name = "integrated_evidence"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_power_tier_thresholds() -> None:
    mod = _load_script()
    assert mod.power_tier(0.9) == "high"
    assert mod.power_tier(0.8) == "high"
    assert mod.power_tier(0.5) == "moderate"
    assert mod.power_tier(0.49) == "low"
    assert mod.power_tier(None) == "unavailable"


def test_assign_labels_definitions() -> None:
    mod = _load_script()
    frame = pl.DataFrame(
        {
            "rsid": ["rs1", "rs2", "rs3", "rs4"],
            "panel_alt": ["A", "C", "G", "T"],
            "gwas_n_traits_q_lt_0_05": [1, 0, 0, 2],
            "eqtl_n_datasets_sig": [0, 3, 0, 1],
            "tsi_differs_fdr05": [False, False, True, True],
        }
    )
    labeled = mod.assign_labels(frame)
    by_rsid = {row["rsid"]: row for row in labeled.iter_rows(named=True)}
    assert by_rsid["rs1"]["evidence_labels"] == mod.LABEL_GWAS
    assert by_rsid["rs2"]["evidence_labels"] == mod.LABEL_EQTL
    assert by_rsid["rs3"]["evidence_labels"] == mod.LABEL_POP
    assert by_rsid["rs4"]["n_evidence_labels"] == 3
    assert mod.LABEL_GWAS in by_rsid["rs4"]["evidence_labels"]
    assert mod.LABEL_EQTL in by_rsid["rs4"]["evidence_labels"]
    assert mod.LABEL_POP in by_rsid["rs4"]["evidence_labels"]


def test_eqtl_direction_consistent_vs_mixed() -> None:
    mod = _load_script()
    panel = pl.DataFrame(
        {
            "rsid": ["rsA", "rsB"],
            "panel_alt": ["A", "C"],
            "gene_symbol": ["G1", "G2"],
            "tsi_differs_fdr05": [False, False],
            "q_tsi_vs_nfe": [1.0, 1.0],
            "testability_rank": [1, 2],
            "power_n1000_or1_5_alpha_bonferroni": [0.9, 0.9],
        }
    )
    assoc = pl.DataFrame(
        {
            "canonical_rsid": ["rsA", "rsA", "rsB", "rsB"],
            "eqtl_alt_allele": ["A", "A", "C", "C"],
            "dataset_id": ["Q1", "Q2", "Q1", "Q2"],
            "beta": [0.2, 0.3, 0.4, -0.5],
            "fdr_bh": [0.01, 0.02, 0.01, 0.01],
        }
    )
    summary = mod.summarize_eqtl(assoc, panel)
    by_rsid = {row["rsid"]: row for row in summary.iter_rows(named=True)}
    assert by_rsid["rsA"]["eqtl_n_datasets_sig"] == 2
    assert by_rsid["rsA"]["eqtl_alt_effect_direction"] == "consistent"
    assert by_rsid["rsB"]["eqtl_n_datasets_sig"] == 2
    assert by_rsid["rsB"]["eqtl_alt_effect_direction"] == "mixed"
