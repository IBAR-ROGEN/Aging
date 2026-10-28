#!/usr/bin/env python3
"""Check that gnomAD NFE and 1000 Genomes CEU refer to the same alt allele.

Flags a row when the absolute NFE−CEU difference is greater than 0.15.
The complement check (NFE closer to 1−CEU than to CEU) is skipped when NFE
AF is between 0.4 and 0.6. In that band the frequency and its complement are
both near 0.5, so a smaller distance to 1−CEU does not show that the alleles
were swapped.

Writes ``analysis/panel/ro_longevity_panel_sanity.tsv`` and records the flag
count on ``analysis/panel/ro_longevity_panel_provenance.json``.

Example:
    uv run python scripts/panel/ro_longevity_panel_sanity.py
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import polars as pl
import typer

from rogen_aging.config import find_repo_root

REPO_ROOT = find_repo_root()
DEFAULT_PANEL = REPO_ROOT / "analysis" / "panel" / "ro_longevity_panel.csv"
DEFAULT_OUT = REPO_ROOT / "analysis" / "panel" / "ro_longevity_panel_sanity.tsv"
DEFAULT_PROVENANCE = REPO_ROOT / "analysis" / "panel" / "ro_longevity_panel_provenance.json"
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


def _update_provenance(path: Path, *, flagged: int, out_path: Path, panel_path: Path) -> None:
    """Add the sanity flag count to the panel provenance JSON."""
    if path.is_file():
        payload: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    else:
        payload = {}
    payload["sanity_check"] = {
        "panel_path": str(panel_path.relative_to(REPO_ROOT))
        if panel_path.is_relative_to(REPO_ROOT)
        else str(panel_path),
        "output_path": str(out_path.relative_to(REPO_ROOT))
        if out_path.is_relative_to(REPO_ROOT)
        else str(out_path),
        "flag_count": flagged,
        "af_gap": AF_GAP,
        "complement_skip_lo": COMPLEMENT_SKIP_LO,
        "complement_skip_hi": COMPLEMENT_SKIP_HI,
        "note": COMPLEMENT_SKIP_NOTE,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main(
    panel: Path = typer.Option(
        DEFAULT_PANEL,
        "--panel",
        help="Panel CSV. Read only.",
    ),
    out: Path = typer.Option(
        DEFAULT_OUT,
        "--out",
        help="Sanity TSV to write (rsid, alt, NFE, CEU, GBR, flags).",
    ),
    provenance: Path = typer.Option(
        DEFAULT_PROVENANCE,
        "--provenance",
        help="Panel provenance JSON to update with the flag count.",
    ),
) -> None:
    """Print NFE, CEU, and GBR alt-allele frequencies and any orientation flags."""
    frame = pl.read_csv(panel)
    typer.echo(COMPLEMENT_SKIP_NOTE)
    flagged = 0
    rows_out: list[dict[str, Any]] = []
    for row in frame.iter_rows(named=True):
        nfe_af = row["af_alt_gnomad_nfe"]
        ceu_af = row["af_alt_1kg_ceu"]
        gbr_af = row["af_alt_1kg_gbr"]
        if nfe_af is None or ceu_af is None:
            typer.echo(f"{row['rsid']} missing NFE or CEU frequency")
            rows_out.append(
                {
                    "rsid": row["rsid"],
                    "alt": row["alt"],
                    "af_alt_gnomad_nfe": nfe_af,
                    "af_alt_1kg_ceu": ceu_af,
                    "af_alt_1kg_gbr": gbr_af,
                    "flags": "missing_nfe_or_ceu",
                    "complement_skipped": False,
                }
            )
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
        rows_out.append(
            {
                "rsid": row["rsid"],
                "alt": row["alt"],
                "af_alt_gnomad_nfe": float(nfe_af),
                "af_alt_1kg_ceu": float(ceu_af),
                "af_alt_1kg_gbr": None if gbr_af is None else float(gbr_af),
                "flags": ",".join(flags) if flags else "",
                "complement_skipped": skipped,
            }
        )
    out.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(rows_out).write_csv(out, separator="\t")
    _update_provenance(provenance, flagged=flagged, out_path=out, panel_path=panel)
    typer.echo(f"flagged {flagged}")
    typer.echo(f"Wrote {out}")
    typer.echo(f"Updated {provenance} sanity_check.flag_count={flagged}")


if __name__ == "__main__":
    typer.run(main)
