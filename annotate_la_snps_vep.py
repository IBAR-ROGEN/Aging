#!/usr/bin/env python3
"""Annotate longevity-associated SNPs (LA-SNPs) via the Ensembl VEP REST API (GRCh38).

Builds a one-row-per-rsID table of inferred functional consequences for manuscript
supplementary tables. Uses per-rsID API calls with local JSON caching and polite
request pacing — no local VCF or bulk download required.
"""

from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Any
from urllib.parse import quote

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# CONSTANTS — edit these paths and tuning knobs before running
# ---------------------------------------------------------------------------

# Curated ~70 LA-SNP set (58 unique rsIDs) from AlphaGenome analysis output.
INPUT_RSID_FILE = Path("alphagenome_impact_analysis.csv")
RSID_COLUMN = "snp"  # column name when INPUT_RSID_FILE is CSV; ignored for plain .txt

OUTPUT_DIR = Path("analysis/vep_annotation")
REST_BASE_URL = "https://rest.ensembl.org"
REQUEST_DELAY_SEC = 0.34
CACHE_DIR = Path("analysis/vep_cache")  # one {rsid}.json file per cached response

REQUEST_TIMEOUT_SEC = 30.0
MAX_RETRIES = 4

# Ensembl VEP impact rank (lower = more severe) for picking SIFT/PolyPhen source transcript.
IMPACT_RANK: dict[str, int] = {
    "HIGH": 0,
    "MODERATE": 1,
    "LOW": 2,
    "MODIFIER": 3,
}

RSID_PATTERN = re.compile(r"^rs\d+$", re.IGNORECASE)


def load_rsids(path: Path, column: str) -> list[str]:
    """Load unique rsIDs from a CSV (column) or plain-text file (one rsID per line)."""
    if not path.is_file():
        raise FileNotFoundError(f"Input rsID file not found: {path.resolve()}")

    suffix = path.suffix.lower()
    if suffix in {".csv", ".tsv"}:
        sep = "\t" if suffix == ".tsv" else ","
        frame = pd.read_csv(path, sep=sep)
        if column not in frame.columns:
            candidates = [column, "snp", "SNP_rsID", "rsID", "rsid"]
            col = next((c for c in candidates if c in frame.columns), None)
            if col is None:
                raise ValueError(
                    f"Column {column!r} not in {path.name}; available: {list(frame.columns)}"
                )
        else:
            col = column
        raw_values = frame[col].tolist()
    else:
        raw_values = path.read_text(encoding="utf-8").splitlines()

    seen: set[str] = set()
    rsids: list[str] = []
    for raw in raw_values:
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            continue
        token = str(raw).strip()
        if not token or not RSID_PATTERN.match(token):
            continue
        key = token.lower()
        if key not in seen:
            seen.add(key)
            rsids.append(token)
    return rsids


def cache_path_for(rsid: str) -> Path:
    """Return the on-disk cache file for a given rsID."""
    safe = quote(rsid, safe="")
    return CACHE_DIR / f"{safe}.json"


def fetch_vep_json(
    rsid: str,
    session: requests.Session,
    *,
    min_interval_sec: float,
    timeout_sec: float,
    max_retries: int,
) -> tuple[list[dict[str, Any]] | None, bool]:
    """Query VEP for one rsID, using cache when present.

    Returns:
        (payload, from_cache) where payload is the parsed JSON list on success,
        or None when the variant is not found / the request ultimately fails.
    """
    cache_file = cache_path_for(rsid)
    if cache_file.is_file():
        with cache_file.open(encoding="utf-8") as handle:
            return json.load(handle), True

    url = f"{REST_BASE_URL}/vep/human/id/{quote(rsid, safe='')}?content-type=application/json"
    attempt = 0
    while True:
        attempt += 1
        try:
            response = session.get(url, timeout=timeout_sec)
        except requests.RequestException as exc:
            print(f"  request error for {rsid}: {exc}")
            return None, False

        if response.status_code == 404:
            return None, False

        if response.status_code in {429, 503}:
            retry_after = response.headers.get("Retry-After")
            if retry_after is not None:
                try:
                    sleep_s = float(retry_after)
                except ValueError:
                    sleep_s = min_interval_sec * (2 ** (attempt - 1))
            else:
                sleep_s = min_interval_sec * (2 ** (attempt - 1))
            sleep_s = max(sleep_s, min_interval_sec)
            if attempt > max_retries:
                print(f"  giving up on {rsid} after HTTP {response.status_code}")
                return None, False
            print(f"  HTTP {response.status_code} for {rsid}; retry in {sleep_s:.2f}s")
            time.sleep(sleep_s)
            continue

        if not response.ok:
            print(f"  HTTP {response.status_code} for {rsid}: {response.text[:200]}")
            return None, False

        payload: list[dict[str, Any]] = response.json()
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        with cache_file.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
        return payload, False


def parse_allele_string(allele_string: str) -> tuple[str, str]:
    """Split Ensembl allele_string (e.g. ``T/A/C``) into ref and alt alleles."""
    alleles = [a.strip() for a in str(allele_string).split("/") if a.strip()]
    if not alleles:
        return "", ""
    ref = alleles[0]
    alts = alleles[1:]
    return ref, ",".join(alts)


def pick_sift_polyphen(
    transcript_consequences: list[dict[str, Any]],
) -> tuple[str, str]:
    """Return SIFT and PolyPhen predictions from the highest-impact transcript that has them."""
    ranked = sorted(
        transcript_consequences,
        key=lambda tc: IMPACT_RANK.get(str(tc.get("impact", "MODIFIER")), 99),
    )
    for tc in ranked:
        sift = tc.get("sift_prediction") or tc.get("sift")
        polyphen = tc.get("polyphen_prediction") or tc.get("polyphen")
        if sift or polyphen:
            return str(sift or ""), str(polyphen or "")
    return "", ""


def extract_annotation(rsid: str, payload: list[dict[str, Any]]) -> dict[str, str] | None:
    """Map a VEP JSON response to one output-table row."""
    if not payload:
        return None

    variant = payload[0]
    transcripts = variant.get("transcript_consequences") or []
    if not isinstance(transcripts, list):
        transcripts = []

    gene_symbols = sorted(
        {
            str(tc.get("gene_symbol")).strip()
            for tc in transcripts
            if tc.get("gene_symbol")
        }
    )
    ref, alt = parse_allele_string(str(variant.get("allele_string", "")))
    sift, polyphen = pick_sift_polyphen(transcripts)

    return {
        "rsID": rsid,
        "chromosome": str(variant.get("seq_region_name", "")),
        "position_GRCh38": str(variant.get("start", "")),
        "ref_allele": ref,
        "alt_allele": alt,
        "most_severe_consequence": str(variant.get("most_severe_consequence", "")),
        "gene_symbols": ";".join(gene_symbols),
        "SIFT": sift,
        "PolyPhen": polyphen,
    }


def build_output_table(rsids: list[str], session: requests.Session) -> tuple[pd.DataFrame, list[str]]:
    """Query (or load from cache) every rsID and assemble the annotation table."""
    rows: list[dict[str, str]] = []
    not_found: list[str] = []
    last_request_end: float | None = None

    for index, rsid in enumerate(rsids, start=1):
        # Pace live API calls; cached hits do not need a delay.
        if last_request_end is not None:
            elapsed = time.monotonic() - last_request_end
            wait = REQUEST_DELAY_SEC - elapsed
            if wait > 0:
                time.sleep(wait)

        print(f"[{index}/{len(rsids)}] {rsid}")
        payload, from_cache = fetch_vep_json(
            rsid,
            session,
            min_interval_sec=REQUEST_DELAY_SEC,
            timeout_sec=REQUEST_TIMEOUT_SEC,
            max_retries=MAX_RETRIES,
        )
        if not from_cache:
            last_request_end = time.monotonic()
        elif last_request_end is None:
            last_request_end = time.monotonic()

        if payload is None:
            not_found.append(rsid)
            continue

        row = extract_annotation(rsid, payload)
        if row is None:
            not_found.append(rsid)
            continue
        rows.append(row)

    return pd.DataFrame(rows), not_found


def main() -> None:
    """Run VEP annotation for all LA-SNPs and write CSV + Excel outputs."""
    rsids = load_rsids(INPUT_RSID_FILE, RSID_COLUMN)
    if not rsids:
        raise SystemExit(f"No rsIDs loaded from {INPUT_RSID_FILE}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_DIR / "la_snp_vep_annotations.csv"
    xlsx_path = OUTPUT_DIR / "la_snp_vep_annotations.xlsx"
    not_found_path = OUTPUT_DIR / "la_snp_vep_not_found.txt"

    session = requests.Session()
    session.headers.update(
        {
            "Content-Type": "application/json",
            "User-Agent": "rogen-aging-vep-annotate/1.0 (ROGEN; academic research; Ensembl REST)",
        }
    )

    table, not_found = build_output_table(rsids, session)
    session.close()

    table.to_csv(csv_path, index=False)
    table.to_excel(xlsx_path, index=False, engine="openpyxl")

    if not_found:
        not_found_path.write_text("\n".join(not_found) + "\n", encoding="utf-8")
        print(f"\nNot found ({len(not_found)}): {', '.join(not_found)}")
    elif not_found_path.is_file():
        not_found_path.unlink()

    n_queried = len(rsids)
    n_annotated = len(table)
    n_missing = len(not_found)

    print("\n--- Summary ---")
    print(f"rsIDs queried:   {n_queried}")
    print(f"rsIDs annotated: {n_annotated}")
    print(f"rsIDs not found: {n_missing}")
    print(f"CSV output:      {csv_path.resolve()}")
    print(f"Excel output:    {xlsx_path.resolve()}")
    if not_found:
        print(f"Not-found list:  {not_found_path.resolve()}")
    print(f"VEP cache dir:   {CACHE_DIR.resolve()}")


if __name__ == "__main__":
    main()
