"""Unit tests for integrated evidence labels, eQTL orientation, and power tiers."""

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
    # Always reload so edits in the same pytest session are picked up.
    sys.modules.pop(name, None)
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
            "eqtl_input_status": ["ok", "ok", "ok", "ok"],
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


def test_missing_eqtl_input_writes_not_available() -> None:
    """Absent associations/annotation must not become 0 or false."""
    mod = _load_script()
    panel = pl.DataFrame(
        {
            "rsid": ["rsMiss"],
            "panel_ref": ["A"],
            "panel_alt": ["G"],
            "gene_symbol": ["GENE"],
            "tsi_differs_fdr05": [False],
            "q_tsi_vs_nfe": [1.0],
            "p_tsi_vs_nfe": [1.0],
            "testability_rank": [1],
            "af_alt_gnomad_nfe": [0.2],
            "af_alt_1kg_tsi": [0.2],
            "power_n1000_or1_5_alpha_bonferroni": [0.9],
        }
    )
    gwas = pl.DataFrame(
        {
            "rsid": ["rsMiss"],
            "panel_alt": ["G"],
            "trait": ["parental lifespan"],
            "direction": ["neutral"],
            "q_bh": [0.5],
            "status": ["present"],
        }
    )
    grid = pl.DataFrame(
        {
            "rsid": ["rsMiss"],
            "n_cases_per_arm": [1000],
            "allelic_or": [1.5],
            "alpha_name": ["bonferroni_0.05_over_47"],
            "power": [0.9],
        }
    )
    table = mod.build_table(
        panel=panel,
        gwas=gwas,
        assoc=None,
        grid=grid,
        eqtl_ok=False,
    )
    row = table.row(0, named=True)
    assert row["eqtl_input_status"] == mod.NOT_AVAILABLE
    assert row["eqtl_n_datasets_sig"] == mod.NOT_AVAILABLE
    assert row["eqtl_alt_effect_direction"] == mod.NOT_AVAILABLE
    assert row["eqtl_datasets_sig"] == mod.NOT_AVAILABLE
    assert row["label_eqtl_support"] == mod.NOT_AVAILABLE
    assert row["eqtl_n_datasets_sig"] != 0
    assert row["label_eqtl_support"] is not False
    assert mod.LABEL_EQTL not in (row["evidence_labels"] or "")


def test_swapped_eqtl_allele_flips_beta_sign() -> None:
    """eQTL effect allele equal to panel ref flips the beta onto panel alt."""
    mod = _load_script()
    panel = pl.DataFrame(
        {
            "rsid": ["rsSwap"],
            "panel_ref": ["A"],
            "panel_alt": ["G"],
            "gene_symbol": ["GENE"],
            "tsi_differs_fdr05": [False],
            "q_tsi_vs_nfe": [1.0],
            "p_tsi_vs_nfe": [1.0],
            "testability_rank": [1],
            "af_alt_gnomad_nfe": [0.25],
            "af_alt_1kg_tsi": [0.25],
            "power_n1000_or1_5_alpha_bonferroni": [0.9],
        }
    )
    # Catalogue reports effect on A (= panel ref); beta +0.4 → panel-alt beta -0.4.
    assoc = pl.DataFrame(
        {
            "canonical_rsid": ["rsSwap"],
            "eqtl_alt_allele": ["A"],
            "dataset_id": ["QTD000146"],
            "beta": [0.4],
            "fdr_bh": [0.01],
            "maf": [0.25],
        }
    )
    summary = mod.summarize_eqtl(assoc, panel)
    row = summary.row(0, named=True)
    assert row["eqtl_n_datasets_sig"] == 1
    assert row["eqtl_n_flipped"] == 1
    assert row["eqtl_alt_effect_direction"] == "consistent"
    # Orient helper alone.
    beta, status, flipped = mod.orient_eqtl_beta(
        panel_ref="A",
        panel_alt="G",
        eqtl_alt="A",
        beta=0.4,
        maf=0.25,
        panel_maf=0.25,
    )
    assert flipped is True
    assert status == "swapped"
    assert beta == -0.4


def test_eqtl_direction_consistent_vs_mixed() -> None:
    mod = _load_script()
    panel = pl.DataFrame(
        {
            "rsid": ["rsA", "rsB"],
            "panel_ref": ["T", "G"],
            "panel_alt": ["A", "C"],
            "gene_symbol": ["G1", "G2"],
            "tsi_differs_fdr05": [False, False],
            "q_tsi_vs_nfe": [1.0, 1.0],
            "p_tsi_vs_nfe": [1.0, 1.0],
            "testability_rank": [1, 2],
            "af_alt_gnomad_nfe": [0.2, 0.2],
            "af_alt_1kg_tsi": [0.2, 0.2],
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
            "maf": [0.2, 0.2, 0.2, 0.2],
        }
    )
    summary = mod.summarize_eqtl(assoc, panel)
    by_rsid = {row["rsid"]: row for row in summary.iter_rows(named=True)}
    assert by_rsid["rsA"]["eqtl_n_datasets_sig"] == 2
    assert by_rsid["rsA"]["eqtl_alt_effect_direction"] == "consistent"
    assert by_rsid["rsB"]["eqtl_n_datasets_sig"] == 2
    assert by_rsid["rsB"]["eqtl_alt_effect_direction"] == "mixed"


def test_palindromic_high_maf_is_ambiguous() -> None:
    mod = _load_script()
    beta, status, flipped = mod.orient_eqtl_beta(
        panel_ref="A",
        panel_alt="T",
        eqtl_alt="A",
        beta=0.3,
        maf=0.45,
        panel_maf=0.45,
    )
    assert status == "ambiguous_palindromic"
    assert beta is None
    assert flipped is False
