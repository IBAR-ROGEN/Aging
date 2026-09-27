"""Detectability calculations for the longevity-variant panel script."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SCRIPT = _REPO_ROOT / "scripts" / "panel" / "ro_longevity_panel_power.py"


def _load_script() -> ModuleType:
    name = "ro_longevity_panel_power"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


panel = _load_script()


def test_missing_variant_file_stops(tmp_path: Path) -> None:
    with pytest.raises(panel.typer.Exit):
        panel.load_variants(tmp_path / "variants_47_input.csv")


def test_wrong_row_count_stops(tmp_path: Path) -> None:
    path = tmp_path / "variants_47_input.csv"
    path.write_text("chrom,pos,ref,alt,rsid,gene_symbol\n1,10,A,G,rs1,GENE\n", encoding="utf-8")
    with pytest.raises(panel.typer.Exit):
        panel.load_variants(path)


def test_maf_is_the_smaller_allele_frequency() -> None:
    assert panel.maf_from_af(0.2) == pytest.approx(0.2)
    assert panel.maf_from_af(0.8) == pytest.approx(0.2)
    assert panel.maf_from_af(None) is None
    assert panel.maf_from_af(1.2) is None


def test_null_odds_ratio_has_power_near_alpha() -> None:
    alpha = 0.01
    power = panel.allelic_power(0.2, 1.0, 1000, alpha)
    assert power is not None
    assert power == pytest.approx(alpha, abs=1e-6)


def test_power_rises_with_sample_size_and_odds_ratio() -> None:
    small = panel.allelic_power(0.2, 1.2, 250, panel.ALPHA_BONFERRONI)
    large_n = panel.allelic_power(0.2, 1.2, 2500, panel.ALPHA_BONFERRONI)
    large_or = panel.allelic_power(0.2, 2.0, 250, panel.ALPHA_BONFERRONI)
    genomewide = panel.allelic_power(0.2, 1.2, 2500, panel.ALPHA_GENOMEWIDE)
    assert small is not None and large_n is not None and large_or is not None
    assert genomewide is not None
    assert large_n > small
    assert large_or > small
    assert genomewide < large_n


def test_missing_or_boundary_maf_has_no_power() -> None:
    assert panel.allelic_power(None, 1.5, 1000, 0.05) is None
    assert panel.allelic_power(0.0, 1.5, 1000, 0.05) is None
    assert panel.allelic_power(1.0, 1.5, 1000, 0.05) is None


def test_ensembl_uses_alt_and_does_not_borrow_another_population() -> None:
    populations = [
        {"population": "1000GENOMES:phase_3:TSI", "allele": "A", "frequency": 0.1},
        {"population": "1000GENOMES:phase_3:TSI", "allele": "G", "frequency": 0.9},
        {"population": "1000GENOMES:phase_3:FIN", "allele": "A", "frequency": 0.4},
        {"population": "gnomAD:NFE", "allele": "A", "frequency": 0.99},
    ]
    frequency, note = panel.alt_allele_frequency(populations, "TSI", "G", "A")
    assert frequency == pytest.approx(0.1)
    assert note is None
    missing, missing_note = panel.alt_allele_frequency(populations, "IBS", "G", "A")
    assert missing is None
    assert missing_note is not None
    complement, complement_note = panel.alt_allele_frequency(
        [{"population": "1000GENOMES:phase_3:CEU", "allele": "G", "frequency": 0.75}],
        "CEU",
        "G",
        "A",
    )
    assert complement == pytest.approx(0.25)
    assert complement_note is None


def test_info_data_release_parser() -> None:
    assert panel.release_number_from_info({"releases": [116]}) == 116
    assert panel.release_number_from_info({"release": 116}) == 116
    assert panel.release_number_from_info({"releases": [115, 116]}) is None
    assert panel.release_number_from_info({"releases": []}) is None


def test_unmatched_alleles_stay_empty() -> None:
    frequency, note = panel.alt_allele_frequency(
        [{"population": "1000GENOMES:phase_3:TSI", "allele": "C", "frequency": 0.2}],
        "TSI",
        "A",
        "G",
    )
    assert frequency is None
    assert note is not None
