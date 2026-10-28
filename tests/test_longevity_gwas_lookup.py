"""Unit tests for longevity GWAS allele harmonization and direction helper."""

from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import polars as pl
import pytest

from rogen_aging.longevity.allele_harmonize import (
    allele_complement,
    harmonize_to_panel_alt,
    is_palindromic,
)

_REPO = Path(__file__).resolve().parent.parent
_SCRIPT = _REPO / "scripts" / "longevity" / "longevity_gwas_lookup.py"


def _load_script() -> ModuleType:
    name = "longevity_gwas_lookup"
    if name in sys.modules:
        del sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_complement_and_palindrome() -> None:
    assert allele_complement("A") == "T"
    assert allele_complement("cg") == "GC"
    assert is_palindromic("A", "T") is True
    assert is_palindromic("G", "A") is False


def test_harmonize_same_swapped_strand() -> None:
    same = harmonize_to_panel_alt(
        panel_ref="G",
        panel_alt="A",
        effect_allele="A",
        other_allele="G",
        beta=0.2,
        eaf=0.3,
        panel_maf=0.3,
    )
    assert same.status == "same"
    assert same.beta_panel_alt == pytest.approx(0.2)

    swapped = harmonize_to_panel_alt(
        panel_ref="G",
        panel_alt="A",
        effect_allele="G",
        other_allele="A",
        beta=0.2,
        eaf=0.3,
        panel_maf=0.3,
    )
    assert swapped.status == "swapped"
    assert swapped.beta_panel_alt == pytest.approx(-0.2)
    assert swapped.eaf_panel_alt == pytest.approx(0.7)

    strand = harmonize_to_panel_alt(
        panel_ref="G",
        panel_alt="A",
        effect_allele="T",
        other_allele="C",
        beta=0.15,
        eaf=0.4,
        panel_maf=0.3,
    )
    assert strand.status == "strand_same"
    assert strand.beta_panel_alt == pytest.approx(0.15)


def test_ambiguous_palindrome_not_flipped() -> None:
    result = harmonize_to_panel_alt(
        panel_ref="A",
        panel_alt="T",
        effect_allele="T",
        other_allele="A",
        beta=0.1,
        eaf=0.5,
        panel_maf=0.45,
    )
    assert result.status == "ambiguous_palindromic"
    assert result.beta_panel_alt is None


def test_benjamini_hochberg_and_direction(tmp_path: Path) -> None:
    lookup = _load_script()
    q = lookup.benjamini_hochberg([0.01, 0.04, None, 0.5])
    assert q[0] is not None and q[0] <= q[1]  # type: ignore[operator]
    assert q[2] is None
    label, note = lookup.direction_vs_reported_longevity_allele(
        panel_ref="G",
        panel_alt="A",
        longevity_associated_allele="A",
        beta_panel_alt=0.2,
    )
    assert label == "concordant"
    assert note == ""
    label, note = lookup.direction_vs_reported_longevity_allele(
        panel_ref="G",
        panel_alt="A",
        longevity_associated_allele=None,
        beta_panel_alt=0.2,
    )
    assert label == "needs_curation"
    assert note == ""


def test_empty_curation_needs_curation_even_if_longevitymap_alt_equals_panel_alt() -> None:
    """longevitymap_alt_allele matching panel alt must not imply concordance.

    Direction uses only the curation file. An empty curated allele is
    needs_curation regardless of what allele_evidence would have said.
    """
    lookup = _load_script()
    # Simulate the September failure mode: LongevityMap alt equals panel alt,
    # GWAS beta is positive, but curation is empty.
    longevitymap_alt_allele = "A"
    panel_alt = "A"
    assert longevitymap_alt_allele == panel_alt
    label, _note = lookup.direction_vs_reported_longevity_allele(
        panel_ref="G",
        panel_alt=panel_alt,
        longevity_associated_allele="",  # empty curation row
        beta_panel_alt=0.25,
    )
    assert label == "needs_curation"
    assert label != "concordant"


def test_direction_strand_check_and_mismatch() -> None:
    lookup = _load_script()
    label, note = lookup.direction_vs_reported_longevity_allele(
        panel_ref="G",
        panel_alt="A",
        longevity_associated_allele="T",  # complement of A
        beta_panel_alt=0.2,
    )
    assert label == "concordant"
    assert note == ""
    # Insertion allele cannot be mapped onto a G/A site even after complement.
    label, note = lookup.direction_vs_reported_longevity_allele(
        panel_ref="G",
        panel_alt="A",
        longevity_associated_allele="GA",
        beta_panel_alt=0.2,
    )
    assert label == "needs_curation"
    assert "after strand check" in note


def test_lazy_rsid_filter_on_gzip_fixture(tmp_path: Path) -> None:
    lookup = _load_script()
    path = tmp_path / "small.tsv.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("rsid\ta1\ta0\tbeta1\tse\tp\tfreq1\n")
        handle.write("rs2032563\tA\tG\t0.1\t0.05\t0.04\t0.2\n")
        handle.write("rs9999999\tC\tT\t0.2\t0.05\t0.01\t0.3\n")
    hits = pl.scan_csv(path, separator="\t").filter(pl.col("rsid").is_in(["rs2032563"])).collect()
    assert hits.height == 1
    assert hits["rsid"][0] == "rs2032563"
    assert lookup.effect_direction(0.1) == "longevity_increasing"
