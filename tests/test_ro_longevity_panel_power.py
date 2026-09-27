"""Detectability calculations for the longevity-variant panel script."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SCRIPT = _REPO_ROOT / "scripts" / "panel" / "ro_longevity_panel_power.py"
_SANITY = _REPO_ROOT / "scripts" / "panel" / "ro_longevity_panel_sanity.py"


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


def _load_sanity() -> ModuleType:
    name = "ro_longevity_panel_sanity"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SANITY)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {_SANITY}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


sanity = _load_sanity()


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


def test_benjamini_hochberg_adjusts_from_the_largest_rank() -> None:
    q_values = panel.benjamini_hochberg([0.01, 0.04, 0.03, None])
    assert q_values[3] is None
    assert q_values[0] == pytest.approx(0.03)
    assert q_values[1] == pytest.approx(0.04)
    assert q_values[2] == pytest.approx(0.04)


def test_count_test_uses_chi_square_only_when_expected_counts_are_large() -> None:
    _small_p, small_method = panel.allele_count_test(1, 5, 2, 8)
    _large_p, large_method = panel.allele_count_test(100, 200, 5000, 10000)
    assert small_method == "fisher"
    assert large_method == "chi_square"


def test_tsi_allele_number_uses_ensembl_counts_or_records_the_sample_size_assumption() -> None:
    counted, counted_n, counted_source = panel.tsi_alt_allele_counts(
        [
            {
                "population": "1000GENOMES:phase_3:TSI",
                "allele": "A",
                "allele_count": 50,
                "frequency": 50 / 214,
            },
            {
                "population": "1000GENOMES:phase_3:TSI",
                "allele": "G",
                "allele_count": 164,
                "frequency": 164 / 214,
            },
        ],
        "G",
        "A",
    )
    assert (counted, counted_n, counted_source) == (50, 214, "ensembl_allele_count")
    assumed, assumed_n, assumed_source = panel.tsi_alt_allele_counts(
        [{"population": "1000GENOMES:phase_3:TSI", "allele": "A", "frequency": 0.25}],
        "G",
        "A",
    )
    assert assumed_source == "assumed_2x_sample_size"
    assert assumed_n == 2 * panel.TSI_PHASE3_N_SAMPLES
    assert assumed == round(0.25 * assumed_n)


def test_nfe_counts_come_from_cache_fields_and_do_not_require_a_new_frequency() -> None:
    assert panel.nfe_counts_from_cache({"af_gnomad_nfe": 0.2, "ac_nfe": 20, "an_nfe": 100}) == (
        20,
        100,
    )
    assert panel.nfe_counts_from_cache({"af_gnomad_nfe": 0.2}) == (None, None)
    assert panel.nfe_counts_from_cache(
        {"af_gnomad_nfe": 0.8582767439809507, "ac_nfe": 0, "an_nfe": 68044}
    ) == (None, None)


def test_difference_wording_names_the_public_panels() -> None:
    assert panel.PANEL_DIFFERENCE == (
        "differs between public reference panels (1000G TSI vs gnomAD NFE)"
    )
    assert "Romanian" not in panel.PANEL_DIFFERENCE
    assert "Southeast" not in panel.PANEL_DIFFERENCE


def test_unmatched_alleles_stay_empty() -> None:
    frequency, note = panel.alt_allele_frequency(
        [{"population": "1000GENOMES:phase_3:TSI", "allele": "C", "frequency": 0.2}],
        "TSI",
        "A",
        "G",
    )
    assert frequency is None
    assert note is not None


def test_complement_check_skips_frequencies_near_one_half() -> None:
    assert sanity.complement_test_applies(0.2) is True
    assert sanity.complement_test_applies(0.4) is False
    assert sanity.complement_test_applies(0.52) is False
    assert sanity.complement_test_applies(0.6) is False
    assert sanity.complement_test_applies(0.61) is True
    near_half = sanity.flag_allele_orientation(0.52, 0.40)
    assert near_half == []
    swapped = sanity.flag_allele_orientation(0.20, 0.75)
    assert swapped == ["abs_nfe_ceu_gt_0.15", "closer_to_ceu_complement"]
    assert "near 0.5" in sanity.COMPLEMENT_SKIP_NOTE
