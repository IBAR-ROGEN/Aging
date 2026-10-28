"""Allele harmonization for GWAS summary statistics vs panel alt.

Same-allele orientation keeps the effect sign. Swapped alleles flip beta and
complement the effect-allele frequency. Strand-flipped alleles are complemented
and compared again. Ambiguous palindromes (A/T or C/G) with panel MAF above
``PALINDROME_MAF_MAX`` are not flipped.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

HarmonizeStatus = Literal[
    "same",
    "swapped",
    "strand_same",
    "strand_swapped",
    "ambiguous_palindromic",
    "allele_mismatch",
]

PALINDROMIC_PAIRS: frozenset[frozenset[str]] = frozenset(
    {
        frozenset({"A", "T"}),
        frozenset({"C", "G"}),
    }
)
PALINDROME_MAF_MAX = 0.42

_COMPLEMENT: dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}


def allele_complement(allele: str) -> str:
    """Return the DNA complement of a single-base or multi-base allele string."""
    parts = [base.upper() for base in allele.strip()]
    if any(base not in _COMPLEMENT for base in parts):
        raise ValueError(f"allele has non-ACGT base: {allele!r}")
    return "".join(_COMPLEMENT[base] for base in parts)


def is_palindromic(ref: str, alt: str) -> bool:
    """True when ref/alt are a single-base A/T or C/G pair."""
    ref_u = ref.strip().upper()
    alt_u = alt.strip().upper()
    if len(ref_u) != 1 or len(alt_u) != 1:
        return False
    return frozenset({ref_u, alt_u}) in PALINDROMIC_PAIRS


def maf_from_af(allele_frequency: float | None) -> float | None:
    """Return the minor-allele frequency of a diallelic frequency, or None."""
    if allele_frequency is None:
        return None
    if not math.isfinite(allele_frequency) or allele_frequency < 0.0 or allele_frequency > 1.0:
        return None
    return float(min(allele_frequency, 1.0 - allele_frequency))


@dataclass(frozen=True)
class HarmonizeResult:
    """Panel-alt-oriented effect after allele matching."""

    status: HarmonizeStatus
    beta_panel_alt: float | None
    eaf_panel_alt: float | None
    flipped: bool


def harmonize_to_panel_alt(
    *,
    panel_ref: str,
    panel_alt: str,
    effect_allele: str,
    other_allele: str,
    beta: float,
    eaf: float | None,
    panel_maf: float | None,
    palindrome_maf_max: float = PALINDROME_MAF_MAX,
) -> HarmonizeResult:
    """Orient a GWAS effect onto the panel alternate allele.

    Args:
        panel_ref: Panel reference allele.
        panel_alt: Panel alternate allele (effect target).
        effect_allele: GWAS effect allele.
        other_allele: GWAS other allele.
        beta: GWAS effect on the effect-allele scale (log OR / log protection).
        eaf: Effect-allele frequency, if present.
        panel_maf: Panel MAF used only for the palindrome guard.
        palindrome_maf_max: MAF above which A/T and C/G sites are ambiguous.

    Returns:
        Harmonized beta/EAF for ``panel_alt``, or an ambiguous/mismatch status.
    """
    ref = panel_ref.strip().upper()
    alt = panel_alt.strip().upper()
    ea = effect_allele.strip().upper()
    oa = other_allele.strip().upper()

    if is_palindromic(ref, alt):
        maf = panel_maf
        if maf is None or maf > palindrome_maf_max:
            return HarmonizeResult(
                status="ambiguous_palindromic",
                beta_panel_alt=None,
                eaf_panel_alt=None,
                flipped=False,
            )

    if ea == alt and oa == ref:
        return HarmonizeResult(status="same", beta_panel_alt=beta, eaf_panel_alt=eaf, flipped=False)
    if ea == ref and oa == alt:
        eaf_flip = None if eaf is None else 1.0 - float(eaf)
        return HarmonizeResult(
            status="swapped",
            beta_panel_alt=-float(beta),
            eaf_panel_alt=eaf_flip,
            flipped=True,
        )

    try:
        ea_c = allele_complement(ea)
        oa_c = allele_complement(oa)
    except ValueError:
        return HarmonizeResult(
            status="allele_mismatch",
            beta_panel_alt=None,
            eaf_panel_alt=None,
            flipped=False,
        )

    if ea_c == alt and oa_c == ref:
        return HarmonizeResult(
            status="strand_same",
            beta_panel_alt=beta,
            eaf_panel_alt=eaf,
            flipped=False,
        )
    if ea_c == ref and oa_c == alt:
        eaf_flip = None if eaf is None else 1.0 - float(eaf)
        return HarmonizeResult(
            status="strand_swapped",
            beta_panel_alt=-float(beta),
            eaf_panel_alt=eaf_flip,
            flipped=True,
        )
    return HarmonizeResult(
        status="allele_mismatch",
        beta_panel_alt=None,
        eaf_panel_alt=None,
        flipped=False,
    )
