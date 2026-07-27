#!/usr/bin/env python3
"""Annotate longevity-associated SNPs (LA-SNPs) with GTEx single-tissue eQTL evidence.

Merges the curated AlphaGenome LA-SNP table with GRCh38 coordinates from the VEP
annotation output, resolves GTEx variant IDs via the GTEx Portal API, and queries
significant single-tissue eQTLs in brain and whole-blood tissues. Raw JSON responses
are cached locally so re-runs skip live API calls.

Uses GTEx Portal API v2 only — no bulk download.
"""

from __future__ import annotations

import hashlib
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
INPUT_IMPACT_FILE = Path("analysis/alphagenome/alphagenome_impact_analysis.csv")
IMPACT_RSID_COLUMN = "snp"

# GRCh38 coordinates from the Ensembl VEP annotation pass (one row per rsID).
INPUT_VEP_FILE = Path("analysis/vep_annotation/la_snp_vep_annotations.xlsx")
VEP_RSID_COLUMN = "rsID"

OUTPUT_DIR = Path("analysis/gtex_annotation")
CACHE_DIR = Path("analysis/gtex_cache")

API_BASE_URL = "https://gtexportal.org/api/v2"
DATASET_ID = "gtex_v10"

# Brain tissues + whole blood (neurodegeneration / peripheral immune angle).
TARGET_TISSUES: list[str] = [
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra",
    "Whole_Blood",
]

REQUEST_DELAY_SEC = 0.5
REQUEST_TIMEOUT_SEC = 30.0
MAX_RETRIES = 4
ITEMS_PER_PAGE = 250

RSID_PATTERN = re.compile(r"^rs\d+$", re.IGNORECASE)

EQTL_OUTPUT_COLUMNS: list[str] = [
    "rsID",
    "gtex_variant_id",
    "gene_symbol",
    "tissue",
    "nes",
    "p_value",
]


def normalize_rsid(raw: object) -> str | None:
    """Return a normalized rsID string, or None if the value is not rs-formatted."""
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return None
    token = str(raw).strip()
    if not token or not RSID_PATTERN.match(token):
        return None
    return token


def load_and_merge_snp_table(
    impact_path: Path,
    impact_rsid_column: str,
    vep_path: Path,
    vep_rsid_column: str,
) -> pd.DataFrame:
    """Load LA-SNPs from the impact table and merge GRCh38 coordinates from VEP output.

    The impact file has 70 gene–SNP rows but only 58 unique rsIDs; we keep one row per
    rsID (first occurrence) and attach VEP coordinates for GTEx variant resolution.
    """
    if not impact_path.is_file():
        raise FileNotFoundError(f"Impact SNP file not found: {impact_path.resolve()}")
    if not vep_path.is_file():
        raise FileNotFoundError(f"VEP coordinate file not found: {vep_path.resolve()}")

    impact = pd.read_csv(impact_path)
    if impact_rsid_column not in impact.columns:
        raise ValueError(
            f"Column {impact_rsid_column!r} missing from {impact_path.name}; "
            f"available: {list(impact.columns)}"
        )

    vep_suffix = vep_path.suffix.lower()
    if vep_suffix == ".csv":
        vep = pd.read_csv(vep_path)
    elif vep_suffix in {".xlsx", ".xls"}:
        vep = pd.read_excel(vep_path)
    else:
        raise ValueError(f"Unsupported VEP file type: {vep_path}")

    if vep_rsid_column not in vep.columns:
        raise ValueError(
            f"Column {vep_rsid_column!r} missing from {vep_path.name}; "
            f"available: {list(vep.columns)}"
        )

    impact["rsID"] = impact[impact_rsid_column].map(normalize_rsid)
    impact = impact.dropna(subset=["rsID"])
    impact_unique = impact.drop_duplicates(subset="rsID", keep="first").copy()

    vep = vep.copy()
    vep["rsID"] = vep[vep_rsid_column].map(normalize_rsid)
    vep = vep.dropna(subset=["rsID"]).drop_duplicates(subset="rsID", keep="first")

    impact_rs = set(impact_unique["rsID"].str.lower())
    vep_rs = set(vep["rsID"].str.lower())
    missing_in_vep = sorted(impact_rs - vep_rs)
    extra_in_vep = sorted(vep_rs - impact_rs)
    if missing_in_vep:
        raise ValueError(
            f"{len(missing_in_vep)} impact rsIDs missing from VEP coordinates: "
            f"{', '.join(missing_in_vep)}"
        )
    if extra_in_vep:
        print(
            f"Note: {len(extra_in_vep)} rsIDs in VEP file not in impact table "
            f"(ignored): {', '.join(extra_in_vep)}"
        )

    coord_cols = ["chromosome", "position_GRCh38", "ref_allele", "alt_allele"]
    missing_coord_cols = [col for col in coord_cols if col not in vep.columns]
    if missing_coord_cols:
        raise ValueError(
            f"VEP file missing coordinate columns: {missing_coord_cols}; "
            f"available: {list(vep.columns)}"
        )

    merged = impact_unique.merge(
        vep[["rsID", *coord_cols]],
        on="rsID",
        how="left",
        validate="1:1",
    )
    if merged["chromosome"].isna().any():
        unresolved = merged.loc[merged["chromosome"].isna(), "rsID"].tolist()
        raise ValueError(f"Merged table has SNPs without coordinates: {unresolved}")

    return merged


def cache_path_for(kind: str, key: str) -> Path:
    """Return a filesystem-safe cache path for a variant or eQTL query."""
    safe = quote(key, safe="")
    if len(safe) > 180:
        digest = hashlib.sha256(key.encode("utf-8")).hexdigest()[:24]
        safe = f"{safe[:80]}__{digest}"
    return CACHE_DIR / f"{kind}_{safe}.json"


def gtex_chromosome(chromosome: object) -> str:
    """Format a chromosome value as GTEx expects (e.g. ``chr2``)."""
    chrom = str(chromosome).strip()
    if chrom.lower().startswith("chr"):
        return chrom if chrom.startswith("chr") else f"chr{chrom[3:]}"
    return f"chr{chrom}"


def fetch_gtex_json(
    endpoint: str,
    params: list[tuple[str, str]],
    cache_file: Path,
    session: requests.Session,
    *,
    min_interval_sec: float,
    timeout_sec: float,
    max_retries: int,
) -> tuple[dict[str, Any] | None, bool]:
    """GET a GTEx API endpoint with caching, pacing, and retry on 429/503."""
    if cache_file.is_file():
        with cache_file.open(encoding="utf-8") as handle:
            return json.load(handle), True

    url = f"{API_BASE_URL}/{endpoint.lstrip('/')}"
    attempt = 0
    while True:
        attempt += 1
        try:
            response = session.get(url, params=params, timeout=timeout_sec)
        except requests.RequestException as exc:
            print(f"  request error for {endpoint}: {exc}")
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
                print(f"  giving up on {endpoint} after HTTP {response.status_code}")
                return None, False
            print(f"  HTTP {response.status_code}; retry in {sleep_s:.2f}s")
            time.sleep(sleep_s)
            continue

        if response.status_code == 400:
            print(f"  HTTP 400 for {endpoint}: {response.text[:300]}")
            return None, False

        if not response.ok:
            print(f"  HTTP {response.status_code} for {endpoint}: {response.text[:300]}")
            return None, False

        payload: dict[str, Any] = response.json()
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        with cache_file.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
        return payload, False


def parse_variant_record(record: dict[str, Any]) -> dict[str, Any]:
    """Extract the GTEx variant fields we care about."""
    return {
        "variantId": str(record.get("variantId", "")),
        "snpId": str(record.get("snpId", "")),
        "chromosome": str(record.get("chromosome", "")),
        "pos": int(record.get("pos", 0)),
        "ref": str(record.get("ref", "")),
        "alt": str(record.get("alt", "")),
    }


def resolve_gtex_variant(
    row: pd.Series,
    session: requests.Session,
    *,
    min_interval_sec: float,
    timeout_sec: float,
    max_retries: int,
) -> dict[str, Any] | None:
    """Resolve a GTEx variant ID for one LA-SNP (rsID first, then GRCh38 locus)."""
    rsid = str(row["rsID"])
    cache_file = cache_path_for("variant", rsid.lower())

    # Primary lookup: dbSNP rsID.
    params: list[tuple[str, str]] = [
        ("snpId", rsid),
        ("datasetId", DATASET_ID),
        ("itemsPerPage", "250"),
    ]
    payload, _from_cache = fetch_gtex_json(
        "dataset/variant",
        params,
        cache_file,
        session,
        min_interval_sec=min_interval_sec,
        timeout_sec=timeout_sec,
        max_retries=max_retries,
    )
    if payload and payload.get("data"):
        variant = parse_variant_record(payload["data"][0])
        expected_pos = int(row["position_GRCh38"])
        if variant["pos"] and variant["pos"] != expected_pos:
            print(
                f"  warning: {rsid} GTEx pos {variant['pos']} != VEP pos {expected_pos}"
            )
        return variant

    # Fallback: chromosome + position from merged VEP coordinates.
    chrom = gtex_chromosome(row["chromosome"])
    pos = str(int(row["position_GRCh38"]))
    fallback_cache = cache_path_for("variant_loc", f"{chrom}_{pos}")
    loc_params: list[tuple[str, str]] = [
        ("chromosome", chrom),
        ("pos", pos),
        ("datasetId", DATASET_ID),
        ("itemsPerPage", "250"),
    ]
    loc_payload, _from_cache = fetch_gtex_json(
        "dataset/variant",
        loc_params,
        fallback_cache,
        session,
        min_interval_sec=min_interval_sec,
        timeout_sec=timeout_sec,
        max_retries=max_retries,
    )
    if not loc_payload or not loc_payload.get("data"):
        return None

    for record in loc_payload["data"]:
        variant = parse_variant_record(record)
        if variant["snpId"].lower() == rsid.lower():
            return variant
        if variant["pos"] == int(row["position_GRCh38"]):
            return variant

    return parse_variant_record(loc_payload["data"][0])


def fetch_significant_eqtls(
    variant_id: str,
    session: requests.Session,
    *,
    min_interval_sec: float,
    timeout_sec: float,
    max_retries: int,
) -> list[dict[str, Any]]:
    """Fetch all significant single-tissue eQTL pages for a GTEx variant in target tissues."""
    tissue_key = "_".join(TARGET_TISSUES)
    cache_file = cache_path_for("eqtl", f"{variant_id}__{tissue_key}")

    if cache_file.is_file():
        with cache_file.open(encoding="utf-8") as handle:
            cached = json.load(handle)
        return cached.get("hits", [])

    hits: list[dict[str, Any]] = []
    page = 0
    number_of_pages = 1

    while page < number_of_pages:
        params: list[tuple[str, str]] = [
            ("variantId", variant_id),
            ("datasetId", DATASET_ID),
            ("itemsPerPage", str(ITEMS_PER_PAGE)),
            ("page", str(page)),
        ]
        for tissue in TARGET_TISSUES:
            params.append(("tissueSiteDetailId", tissue))

        page_cache = cache_path_for(
            "eqtl_page",
            f"{variant_id}__{tissue_key}__p{page}",
        )
        payload, _from_cache = fetch_gtex_json(
            "association/singleTissueEqtl",
            params,
            page_cache,
            session,
            min_interval_sec=min_interval_sec,
            timeout_sec=timeout_sec,
            max_retries=max_retries,
        )
        if payload is None:
            break

        data = payload.get("data") or []
        if isinstance(data, list):
            hits.extend(data)

        paging = payload.get("paging_info") or {}
        number_of_pages = int(paging.get("numberOfPages", 1))
        page += 1

    target_set = set(TARGET_TISSUES)
    hits = [hit for hit in hits if hit.get("tissueSiteDetailId") in target_set]

    cache_file.parent.mkdir(parents=True, exist_ok=True)
    with cache_file.open("w", encoding="utf-8") as handle:
        json.dump({"variantId": variant_id, "hits": hits}, handle, indent=2)

    return hits


def eqtl_hit_to_row(rsid: str, gtex_variant_id: str, hit: dict[str, Any]) -> dict[str, Any]:
    """Map one GTEx singleTissueEqtl record to an output-table row."""
    return {
        "rsID": rsid,
        "gtex_variant_id": gtex_variant_id,
        "gene_symbol": str(hit.get("geneSymbol", "")),
        "tissue": str(hit.get("tissueSiteDetailId", "")),
        "nes": float(hit.get("nes", float("nan"))),
        "p_value": float(hit.get("pValue", float("nan"))),
    }


def build_eqtl_table(
    snp_table: pd.DataFrame,
    session: requests.Session,
) -> tuple[pd.DataFrame, list[str], dict[str, str]]:
    """Resolve variants and collect eQTL hits for every LA-SNP."""
    rows: list[dict[str, Any]] = []
    unresolved: list[str] = []
    resolved_variant_ids: dict[str, str] = {}
    last_request_end: float | None = None

    for index, (_, snp_row) in enumerate(snp_table.iterrows(), start=1):
        rsid = str(snp_row["rsID"])
        print(f"[{index}/{len(snp_table)}] {rsid}")

        if last_request_end is not None:
            elapsed = time.monotonic() - last_request_end
            wait = REQUEST_DELAY_SEC - elapsed
            if wait > 0:
                time.sleep(wait)

        variant = resolve_gtex_variant(
            snp_row,
            session,
            min_interval_sec=REQUEST_DELAY_SEC,
            timeout_sec=REQUEST_TIMEOUT_SEC,
            max_retries=MAX_RETRIES,
        )
        last_request_end = time.monotonic()

        if variant is None or not variant.get("variantId"):
            unresolved.append(rsid)
            continue

        gtex_variant_id = variant["variantId"]
        resolved_variant_ids[rsid] = gtex_variant_id
        print(f"  variantId: {gtex_variant_id}")

        if last_request_end is not None:
            elapsed = time.monotonic() - last_request_end
            wait = REQUEST_DELAY_SEC - elapsed
            if wait > 0:
                time.sleep(wait)

        hits = fetch_significant_eqtls(
            gtex_variant_id,
            session,
            min_interval_sec=REQUEST_DELAY_SEC,
            timeout_sec=REQUEST_TIMEOUT_SEC,
            max_retries=MAX_RETRIES,
        )
        last_request_end = time.monotonic()

        if hits:
            print(f"  eQTL hits in target tissues: {len(hits)}")
        else:
            print("  no significant eQTLs in target tissues")

        for hit in hits:
            rows.append(eqtl_hit_to_row(rsid, gtex_variant_id, hit))

    table = pd.DataFrame(rows, columns=EQTL_OUTPUT_COLUMNS)
    return table, unresolved, resolved_variant_ids


def main() -> None:
    """Run GTEx eQTL annotation for all LA-SNPs and write CSV + Excel outputs."""
    snp_table = load_and_merge_snp_table(
        INPUT_IMPACT_FILE,
        IMPACT_RSID_COLUMN,
        INPUT_VEP_FILE,
        VEP_RSID_COLUMN,
    )
    n_snps = len(snp_table)
    print(f"Merged LA-SNP table: {n_snps} unique rsIDs with GRCh38 coordinates")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_DIR / "la_snp_gtex_eqtls.csv"
    xlsx_path = OUTPUT_DIR / "la_snp_gtex_eqtls.xlsx"
    unresolved_path = OUTPUT_DIR / "la_snp_gtex_unresolved.txt"

    session = requests.Session()
    session.headers.update(
        {
            "Accept": "application/json",
            "User-Agent": (
                "rogen-aging-gtex-annotate/1.0 (ROGEN; academic research; GTEx Portal API)"
            ),
        }
    )

    eqtl_table, unresolved, resolved = build_eqtl_table(snp_table, session)
    session.close()

    eqtl_table.to_csv(csv_path, index=False)
    eqtl_table.to_excel(xlsx_path, index=False, engine="openpyxl")

    if unresolved:
        unresolved_path.write_text("\n".join(unresolved) + "\n", encoding="utf-8")
    elif unresolved_path.is_file():
        unresolved_path.unlink()

    snps_with_hits = eqtl_table["rsID"].nunique() if not eqtl_table.empty else 0
    total_hits = len(eqtl_table)

    print("\n--- Summary ---")
    print(f"LA-SNPs queried:              {n_snps}")
    print(f"Variants resolved:            {len(resolved)}")
    print(f"SNPs with >=1 target eQTL:    {snps_with_hits}")
    print(f"Total eQTL hits (long table): {total_hits}")
    print(f"Unresolved SNPs:              {len(unresolved)}")
    if unresolved:
        print(f"  {', '.join(unresolved)}")
    print(f"CSV output:                   {csv_path.resolve()}")
    print(f"Excel output:                 {xlsx_path.resolve()}")
    if unresolved:
        print(f"Unresolved list:              {unresolved_path.resolve()}")
    print(f"GTEx cache dir:               {CACHE_DIR.resolve()}")


if __name__ == "__main__":
    main()
