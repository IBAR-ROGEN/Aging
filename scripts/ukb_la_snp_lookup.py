#!/usr/bin/env python3
"""Offline UK Biobank (UKB) genotype extraction manifest builder for ROGEN.

This module reads a gene–SNP overlap table from Excel, resolves each dbSNP
rs identifier to a GRCh38 chromosomal locus via the Ensembl Variation REST
API, and writes a CSV manifest suitable for downstream bulk genotype
extraction (for example, imputed genotype resources often discussed in
conjunction with UK Biobank showcase field **22418** — genotype/imputed
bulk files; this script does **not** call DNAnexus, dx-toolkit, or dxFUSE).

The manifest encodes chromosome and GRCh38 position for harmonisation with
contemporary reference-based workflows. Native UKB imputed BGEN products
are historically GRCh37-aligned; confirm liftover requirements before
matching coordinates to released UKB variant sets.

Example:
    uv run python scripts/ukb_la_snp_lookup.py \\
        --input overlapping_genes_with_snps.xlsx \\
        --output analysis/ukb_snp_manifest_v0.1.csv
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import quote

import pandas as pd
import requests

ENSEMBL_VARIATION_BASE = "https://rest.ensembl.org/variation/human"
DEFAULT_ASSEMBLY = "GRCh38"
# Conservative pacing; Ensembl recommends modest parallelism and polite usage.
DEFAULT_MIN_INTERVAL_SEC = 0.34
DEFAULT_TIMEOUT_SEC = 30.0
DEFAULT_MAX_RETRIES = 4

LOG = logging.getLogger(__name__)


@dataclass(frozen=True)
class Grch38Locus:
    """Chromosomal coordinates on the GRCh38 primary assembly."""

    chromosome: str
    position: int


def _pick_grch38_chromosome_mapping(mappings: list[dict[str, Any]]) -> dict[str, Any] | None:
    """Select the preferred GRCh38 mapping dict from Ensembl ``mappings`` list.

    Prefers entries on the primary ``chromosome`` coordinate system and
    ``GRCh38`` assembly, excluding patch/scaffold contigs where a primary
    chromosomal placement exists.

    Args:
        mappings: The ``mappings`` array from a Ensembl variation JSON payload.

    Returns:
        The chosen mapping dictionary, or ``None`` if none qualifies.
    """
    chromosomal_grch38: list[dict[str, Any]] = []
    for m in mappings:
        if m.get("assembly_name") != DEFAULT_ASSEMBLY:
            continue
        if m.get("coord_system") != "chromosome":
            continue
        seq_region = str(m.get("seq_region_name", ""))
        if not seq_region or seq_region.upper().endswith("_PATCH"):
            continue
        chromosomal_grch38.append(m)

    if not chromosomal_grch38:
        return None
    # Prefer shortest region name that looks like 1–22, X, Y, MT (stable sort).
    chromosomal_grch38.sort(key=lambda x: (len(str(x.get("seq_region_name", ""))), str(x.get("seq_region_name", ""))))
    return chromosomal_grch38[0]


def _mapping_to_locus(mapping: dict[str, Any]) -> Grch38Locus | None:
    """Convert a single Ensembl mapping dict to :class:`Grch38Locus`."""
    chrom = mapping.get("seq_region_name")
    start = mapping.get("start")
    end = mapping.get("end")
    if chrom is None or start is None:
        return None
    try:
        pos = int(start)
    except (TypeError, ValueError):
        return None
    if end is not None:
        try:
            end_i = int(end)
        except (TypeError, ValueError):
            end_i = pos
        if end_i != pos:
            LOG.warning(
                "Variant spans multiple bases (start=%s end=%s); using start as SNP position.",
                start,
                end,
            )
    return Grch38Locus(chromosome=str(chrom), position=pos)


def query_ensembl_rsids_grch38(
    rs_ids: Iterable[str],
    *,
    session: requests.Session | None = None,
    min_interval_sec: float = DEFAULT_MIN_INTERVAL_SEC,
    timeout_sec: float = DEFAULT_TIMEOUT_SEC,
    max_retries: int = DEFAULT_MAX_RETRIES,
) -> dict[str, Grch38Locus | None]:
    """Resolve rs identifiers to GRCh38 loci using the Ensembl Variation REST API.

    For each distinct rs identifier, performs ``GET /variation/human/{id}``
    (see Ensembl REST documentation: variation by species and identifier).
    Applies sequential requests with a minimum wall-clock interval between
    successful attempts (rate limiting), retries transient HTTP failures and
    ``429 Too Many Requests`` with exponential backoff, and respects
    ``Retry-After`` when provided.

    Args:
        rs_ids: Iterable of dbSNP-style identifiers (e.g. ``rs123``). Empty
            or whitespace-only entries are skipped.
        session: Optional :class:`requests.Session` for connection pooling.
        min_interval_sec: Minimum seconds between completed request attempts
            (post-retry success path included for simplicity).
        timeout_sec: Per-request socket timeout.
        max_retries: Maximum retry attempts per identifier for retryable
            status codes (429, 502, 503, 504).

    Returns:
        Mapping from normalised rs identifier string to :class:`Grch38Locus`
        on success, or ``None`` if the variant is unknown or lacks a usable
        GRCh38 chromosomal mapping.
    """
    own_session = session is None
    sess = session or requests.Session()
    sess.headers.setdefault(
        "Content-Type",
        "application/json",
    )
    sess.headers.setdefault(
        "User-Agent",
        "rogen-aging-ukb-manifest/0.1 (ROGEN; academic research; Ensembl REST)",
    )

    results: dict[str, Grch38Locus | None] = {}
    last_request_end: float | None = None

    def pace() -> None:
        nonlocal last_request_end
        if last_request_end is None:
            return
        elapsed = time.monotonic() - last_request_end
        wait = min_interval_sec - elapsed
        if wait > 0:
            time.sleep(wait)

    unique_norm: list[str] = []
    seen: set[str] = set()
    for raw in rs_ids:
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            continue
        s = str(raw).strip()
        if not s:
            continue
        key = s
        if key not in seen:
            seen.add(key)
            unique_norm.append(key)

    for rs in unique_norm:
        pace()
        url = f"{ENSEMBL_VARIATION_BASE}/{quote(rs, safe='')}"
        attempt = 0
        while True:
            attempt += 1
            try:
                resp = sess.get(url, timeout=timeout_sec)
            except requests.RequestException as exc:
                LOG.error("Request failed for %s: %s", rs, exc)
                results[rs] = None
                last_request_end = time.monotonic()
                break

            if resp.status_code == 404:
                LOG.warning("Ensembl returned 404 for %s (unknown or retired id).", rs)
                results[rs] = None
                last_request_end = time.monotonic()
                break

            if resp.status_code == 429 or resp.status_code in (502, 503, 504):
                retry_after_hdr = resp.headers.get("Retry-After")
                if retry_after_hdr is not None:
                    try:
                        sleep_s = float(retry_after_hdr)
                    except ValueError:
                        sleep_s = min_interval_sec * (2 ** (attempt - 1))
                else:
                    sleep_s = min_interval_sec * (2 ** (attempt - 1))
                sleep_s = max(sleep_s, min_interval_sec)
                if attempt > max_retries:
                    LOG.error(
                        "Giving up on %s after %s attempts (last status=%s).",
                        rs,
                        max_retries,
                        resp.status_code,
                    )
                    results[rs] = None
                    last_request_end = time.monotonic()
                    break
                LOG.warning(
                    "HTTP %s for %s; sleeping %.2fs then retrying (%s/%s).",
                    resp.status_code,
                    rs,
                    sleep_s,
                    attempt,
                    max_retries,
                )
                time.sleep(sleep_s)
                continue

            if not resp.ok:
                LOG.error(
                    "HTTP %s for %s: %s",
                    resp.status_code,
                    rs,
                    resp.text[:500],
                )
                results[rs] = None
                last_request_end = time.monotonic()
                break

            try:
                payload: dict[str, Any] = resp.json()
            except json.JSONDecodeError as exc:
                LOG.error("Invalid JSON for %s: %s", rs, exc)
                results[rs] = None
                last_request_end = time.monotonic()
                break

            mappings = payload.get("mappings")
            if not isinstance(mappings, list):
                LOG.warning("No mappings list in Ensembl payload for %s.", rs)
                results[rs] = None
                last_request_end = time.monotonic()
                break

            chosen = _pick_grch38_chromosome_mapping(mappings)
            if chosen is None:
                LOG.warning("No GRCh38 chromosomal mapping for %s.", rs)
                results[rs] = None
                last_request_end = time.monotonic()
                break

            locus = _mapping_to_locus(chosen)
            results[rs] = locus
            last_request_end = time.monotonic()
            break

    if own_session:
        sess.close()

    return results


def ukb_expected_chunk(chromosome: str, position: int) -> str:
    """Derive a human-readable UKB-oriented chunk label for imputed bulk planning.

    UK Biobank distributes imputed genotype data in large chromosome-scoped
    bulk files (resource families often referenced alongside showcase field
    **22418**). This function does **not** compute DNAnexus file handles; it
    records a deterministic string that teams can use to batch extractions
    by chromosome while retaining GRCh38 coordinates from Ensembl.

    Args:
        chromosome: Sequence region name from Ensembl (e.g. ``9``, ``X``).
        position: One-based GRCh38 position on that chromosome.

    Returns:
        A short chunk descriptor tying field exemplar, chromosome, and locus.
    """
    return f"F22418_imputed_chr{chromosome}_GRCh38_{position}"


def read_overlap_table(path: Path) -> pd.DataFrame:
    """Load the gene–SNP overlap spreadsheet.

    Args:
        path: Path to ``overlapping_genes_with_snps.xlsx`` (or compatible).

    Returns:
        DataFrame with at least ``Gene`` and ``SNP_rsID`` columns.

    Raises:
        ValueError: If required columns are missing.
    """
    df = pd.read_excel(path, engine="openpyxl")
    required = {"Gene", "SNP_rsID"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input Excel is missing required columns: {sorted(missing)}")
    return df


def build_manifest(
    df: pd.DataFrame,
    *,
    min_interval_sec: float,
    timeout_sec: float,
    max_retries: int,
) -> pd.DataFrame:
    """Attach GRCh38 coordinates and UKB chunk labels to overlap rows.

    Args:
        df: Output of :func:`read_overlap_table`.
        min_interval_sec: Passed to :func:`query_ensembl_rsids_grch38`.
        timeout_sec: Passed to :func:`query_ensembl_rsids_grch38`.
        max_retries: Passed to :func:`query_ensembl_rsids_grch38`.

    Returns:
        DataFrame with columns ``Gene``, ``SNP_rsID``, ``Chromosome``,
        ``Position_GRCh38``, ``UKB_Expected_Chunk``.
    """
    work = df.copy()
    work["Gene"] = work["Gene"].map(lambda x: str(x).strip() if pd.notna(x) else "")
    work["SNP_rsID"] = work["SNP_rsID"].map(lambda x: str(x).strip() if pd.notna(x) else "")
    rs_to_locus = query_ensembl_rsids_grch38(
        work["SNP_rsID"].tolist(),
        min_interval_sec=min_interval_sec,
        timeout_sec=timeout_sec,
        max_retries=max_retries,
    )

    chromosomes: list[str] = []
    positions: list[int | None] = []
    chunks: list[str] = []

    for rs in work["SNP_rsID"]:
        if not rs:
            chromosomes.append("")
            positions.append(None)
            chunks.append("")
            continue
        loc = rs_to_locus.get(rs)
        if loc is None:
            chromosomes.append("")
            positions.append(None)
            chunks.append("")
            continue
        chromosomes.append(loc.chromosome)
        positions.append(loc.position)
        chunks.append(ukb_expected_chunk(loc.chromosome, loc.position))

    out = pd.DataFrame(
        {
            "Gene": work["Gene"],
            "SNP_rsID": work["SNP_rsID"],
            "Chromosome": chromosomes,
            "Position_GRCh38": positions,
            "UKB_Expected_Chunk": chunks,
        }
    )
    return out


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Build an offline CSV manifest of GRCh38 loci for UKB genotype "
            "extraction planning from a gene–SNP overlap Excel file."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("overlapping_genes_with_snps.xlsx"),
        help="Path to overlapping_genes_with_snps.xlsx",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/ukb_snp_manifest_v0.1.csv"),
        help="Output CSV manifest path",
    )
    parser.add_argument(
        "--min-interval",
        type=float,
        default=DEFAULT_MIN_INTERVAL_SEC,
        help="Minimum seconds between Ensembl HTTP request completions",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=DEFAULT_TIMEOUT_SEC,
        help="Per-request timeout in seconds",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=DEFAULT_MAX_RETRIES,
        help="Maximum retries per rs ID for transient HTTP failures",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging verbosity",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point: read overlap table, query Ensembl, write manifest CSV."""
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, str(args.log_level)),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    input_path: Path = args.input
    output_path: Path = args.output

    if not input_path.is_file():
        LOG.error("Input file does not exist: %s", input_path.resolve())
        return 1

    LOG.info("Reading overlap table from %s", input_path.resolve())
    df = read_overlap_table(input_path)
    LOG.info("Loaded %s rows; resolving distinct rs identifiers via Ensembl", len(df))

    manifest = build_manifest(
        df,
        min_interval_sec=float(args.min_interval),
        timeout_sec=float(args.timeout),
        max_retries=int(args.max_retries),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(output_path, index=False)
    LOG.info("Wrote manifest (%s rows) to %s", len(manifest), output_path.resolve())
    return 0


if __name__ == "__main__":
    sys.exit(main())
