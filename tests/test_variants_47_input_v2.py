"""The versioned 47-variant input must use a gnomAD NFE allele."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_VARIANTS = _REPO_ROOT / "data" / "processed" / "variants_47_input_v2.csv"
_EVIDENCE = _REPO_ROOT / "data" / "processed" / "allele_evidence_47.csv"
_NFE_AF_MIN = 0.001


def _allele_frequencies(field: str) -> dict[str, float]:
    frequencies: dict[str, float] = {}
    if not field.strip():
        return frequencies
    for part in field.split("|"):
        allele, value = part.split("=")
        frequencies[allele] = float(value)
    return frequencies


def test_v2_alts_are_in_gnomad_nfe_and_do_not_contradict_longevitymap() -> None:
    """Fail when an input alt is missing from gnomAD NFE or disagrees with LongevityMap.

    LongevityMap contradicts the input when its ``state = alt`` allele is also
    a gnomAD v4 NFE allele with frequency above 0.001 and the input alt is a
    different base. A LongevityMap allele that is absent or near-absent in NFE
    is not treated as the allele to enforce.
    """
    with _VARIANTS.open(newline="", encoding="utf-8") as handle:
        variants = list(csv.DictReader(handle))
    with _EVIDENCE.open(newline="", encoding="utf-8") as handle:
        evidence = {row["rsid"]: row for row in csv.DictReader(handle)}
    assert len(variants) == 47
    failures: list[str] = []
    for row in variants:
        rsid = row["rsid"]
        alt = row["alt"]
        recorded = evidence[rsid]
        common = _allele_frequencies(recorded["gnomad_nfe_af_gt_0_001"])
        af = common.get(alt)
        if af is None or af <= _NFE_AF_MIN:
            failures.append(f"{rsid} alt {alt} is absent from gnomAD NFE above 0.001")
        longevitymap = {
            allele for allele in recorded["longevitymap_alt_allele"].split("|") if allele
        }
        common_longevitymap = longevitymap & set(common)
        if common_longevitymap and alt not in common_longevitymap:
            failures.append(
                f"{rsid} alt {alt} contradicts LongevityMap {sorted(common_longevitymap)}"
            )
    if failures:
        pytest.fail("; ".join(failures))
