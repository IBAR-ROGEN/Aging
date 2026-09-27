#!/usr/bin/env python3
"""Check that gnomAD NFE and 1000 Genomes CEU refer to the same alt allele.

Flags a row when the absolute NFE−CEU difference is greater than 0.15.
The complement check (NFE closer to 1−CEU than to CEU) is skipped when NFE
AF is between 0.4 and 0.6. In that band the frequency and its complement are
both near 0.5, so a smaller distance to 1−CEU does not show that the alleles
were swapped.

Example:
    uv run python scripts/panel/ro_longevity_panel_sanity.py
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import typer

from rogen_aging.config import find_repo_root

REPO_ROOT = find_repo_root()
DEFAULT_PANEL = REPO_ROOT / "analysis" / "panel" / "ro_longevity_panel.csv"
COMPLEMENT_SKIP_LO = 0.4
COMPLEMENT_SKIP_HI = 0.6
AF_GAP = 0.15
COMPLEMENT_SKIP_NOTE = (
    "Complement test skipped when NFE AF is between 0.4 and 0.6: the allele "
    "frequency and its complement are both near 0.5, so a smaller distance to "
    "1 - CEU does not show that the alleles were swapped."
)


def complement_test_applies(nfe_af: float) -> bool:
    """Return false when NFE AF is too close to 0.5 for a complement check."""
    return not (COMPLEMENT_SKIP_LO <= nfe_af <= COMPLEMENT_SKIP_HI)


def flag_allele_orientation(nfe_af: float, ceu_af: float) -> list[str]:
    """Return orientation flags for one alt-allele pair. Frequencies are not changed."""
    flags: list[str] = []
    if abs(nfe_af - ceu_af) > AF_GAP:
        flags.append("abs_nfe_ceu_gt_0.15")
    if complement_test_applies(nfe_af) and abs(nfe_af - (1.0 - ceu_af)) < abs(nfe_af - ceu_af):
        flags.append("closer_to_ceu_complement")
    return flags


def main(
    panel: Path = typer.Option(
        DEFAULT_PANEL,
        "--panel",
        help="Panel CSV. Read only.",
    ),
) -> None:
    """Print NFE, CEU, and GBR alt-allele frequencies and any orientation flags."""
    frame = pl.read_csv(panel)
    typer.echo(COMPLEMENT_SKIP_NOTE)
    flagged = 0
    for row in frame.iter_rows(named=True):
        nfe_af = row["af_alt_gnomad_nfe"]
        ceu_af = row["af_alt_1kg_ceu"]
        gbr_af = row["af_alt_1kg_gbr"]
        if nfe_af is None or ceu_af is None:
            typer.echo(f"{row['rsid']} missing NFE or CEU frequency")
            continue
        flags = flag_allele_orientation(float(nfe_af), float(ceu_af))
        skipped = not complement_test_applies(float(nfe_af))
        mark = ""
        if flags:
            flagged += 1
            mark = " FLAG " + ",".join(flags)
        elif skipped:
            mark = " complement_skipped"
        gbr_text = "" if gbr_af is None else f"{float(gbr_af):.6f}"
        typer.echo(
            f"{row['rsid']} alt={row['alt']} NFE={float(nfe_af):.6f} "
            f"CEU={float(ceu_af):.6f} GBR={gbr_text}{mark}"
        )
    typer.echo(f"flagged {flagged}")


if __name__ == "__main__":
    typer.run(main)
