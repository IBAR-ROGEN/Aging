#!/usr/bin/env python3
"""
Search public databases for NGS sequencing datasets relevant to aging.

Uses NCBI SRA (Sequence Read Archive) via E-utilities. Results are annotated
with whether each dataset was generated with Oxford Nanopore or another platform.

Output: CSV table with accession, title, platform, is_oxford_nanopore, and related fields.

Usage:
    python search_aging_ngs_datasets.py [--max-results N] [--output FILE]
    uv run python search_aging_ngs_datasets.py

Requires: requests (in pyproject.toml). Optional: NCBI_API_KEY in .env for higher rate limits.
"""

import argparse
import csv
import os
import sys
import time
from pathlib import Path
from xml.etree import ElementTree as ET

import requests

# NCBI E-utilities base URL
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Search terms for aging-related NGS studies (combined with OR)
AGING_TERMS = [
    "aging",
    "ageing",
    "longevity",
    "senescence",
    "biological age",
    "epigenetic age",
    "DNA methylation age",
]

# Platform strings that indicate Oxford Nanopore (case-insensitive)
OXFORD_NANOPORE_INDICATORS = (
    "oxford nanopore",
    "nanopore",
    "minion",
    "promethion",
    "gridion",
    "flongle",
)


def get_ncbi_api_key() -> str | None:
    """Return NCBI API key from environment (optional; improves rate limits)."""
    return os.environ.get("NCBI_API_KEY") or os.environ.get("NCBI_API_KEY_here")


def is_oxford_nanopore(platform_str: str | None) -> bool:
    """Return True if the platform string indicates Oxford Nanopore sequencing."""
    if not platform_str or not isinstance(platform_str, str):
        return False
    lower = platform_str.lower().strip()
    return any(ind in lower for ind in OXFORD_NANOPORE_INDICATORS)


def normalize_platform(platform_str: str | None) -> str:
    """Return a short, consistent platform label."""
    if not platform_str or not isinstance(platform_str, str):
        return "unknown"
    s = platform_str.strip()
    if is_oxford_nanopore(s):
        return "Oxford Nanopore"
    # Common SRA platform values
    if "illumina" in s.lower():
        return "Illumina"
    if "pacbio" in s.lower() or "pac bio" in s.lower():
        return "PacBio"
    if "ion torrent" in s.lower():
        return "Ion Torrent"
    if "bgiseq" in s.lower() or "mgi" in s.lower():
        return "MGI/BGI"
    return s or "unknown"


def search_sra(
    query: str,
    max_results: int = 100,
    api_key: str | None = None,
) -> list[str]:
    """Return list of SRA UIDs (numeric) for the given query."""
    params = {
        "db": "sra",
        "term": query,
        "retmax": min(max_results, 500),
        "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key

    r = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    data = r.json()
    id_list = data.get("esearchresult", {}).get("idlist", [])
    return id_list


def _sra_link(accession: str, uid: str) -> str:
    """Build NCBI SRA link for the given accession or UID."""
    if accession and str(accession).strip() and not str(accession).isdigit():
        return f"https://www.ncbi.nlm.nih.gov/sra/{accession.strip()}"
    if uid:
        return f"https://www.ncbi.nlm.nih.gov/sra/?term={uid}"
    return ""


def _description_from_items(item_map: dict, title: str, max_len: int = 2000) -> str:
    """Extract description from SRA summary items (Summary, Abstract, Design, Title)."""
    raw = (
        item_map.get("Summary")
        or item_map.get("Abstract")
        or item_map.get("Design")
        or item_map.get("title")
        or title
    )
    if not raw or not isinstance(raw, str):
        return ""
    s = raw.strip()
    return s[:max_len] + ("..." if len(s) > max_len else "")


def _number_of_samples_from_doc(doc: ET.Element, item_map: dict, accession: str) -> str:
    """Infer number of samples from SRA DocSum (RunCount, total_runs, or count of Runs list)."""
    # Explicit count fields
    run_count = item_map.get("RunCount") or item_map.get("total_runs") or item_map.get("n_runs")
    if run_count is not None and str(run_count).strip().isdigit():
        return str(int(run_count))
    # Count Run items in a List (e.g. Item Name="Runs" Type="List" with Item children)
    for item in doc.findall("Item"):
        if item.get("Type") == "List" and (item.get("Name") or "").lower() in ("runs", "run"):
            children = item.findall("Item")
            if children:
                return str(len(children))
    # Single run/experiment record → 1 sample
    if (item_map.get("Run") or accession or "").strip():
        return "1"
    return ""


def fetch_sra_summaries(
    uid_list: list[str],
    api_key: str | None = None,
) -> list[dict]:
    """Fetch SRA document summaries for the given UIDs. Returns list of dicts with accession, title, platform, etc."""
    if not uid_list:
        return []

    params = {
        "db": "sra",
        "id": ",".join(uid_list),
        "retmode": "xml",
    }
    if api_key:
        params["api_key"] = api_key

    r = requests.get(f"{EUTILS_BASE}/esummary.fcgi", params=params, timeout=60)
    r.raise_for_status()
    root = ET.fromstring(r.content)

    results = []
    for doc in root.findall(".//DocSum"):
        uid_el = doc.find("Id")
        uid = uid_el.text if uid_el is not None else ""
        item_map = {}
        for item in doc.findall("Item"):
            name = item.get("Name")
            if name is not None and item.get("Type") != "List":
                item_map[name] = item.text or ""

        # Prefer Run accession (SRR), then Experiment (SRX), then Study (SRP)
        accession = (
            item_map.get("Run") or item_map.get("Accession") or item_map.get("Experiment") or item_map.get("Study") or uid
        )
        title = item_map.get("Title") or item_map.get("Run") or ""
        platform = item_map.get("Platform") or item_map.get("PlatformInstrument") or ""
        organism = item_map.get("Organism") or item_map.get("OrganismScientificName") or ""
        study = item_map.get("Study") or item_map.get("StudyAcc") or ""
        if not platform and "Platform" in str(item_map):
            platform = item_map.get("Platform", "")

        description = _description_from_items(item_map, title)
        number_of_samples = _number_of_samples_from_doc(doc, item_map, accession)
        link = _sra_link(accession, uid)

        results.append({
            "uid": uid,
            "accession": accession,
            "title": title,
            "description": description,
            "platform_raw": platform,
            "platform": normalize_platform(platform),
            "is_oxford_nanopore": is_oxford_nanopore(platform),
            "organism": organism,
            "study_accession": study,
            "number_of_samples": number_of_samples,
            "link": link,
        })
    return results


def build_search_query(
    extra_terms: list[str] | None = None,
    nanopore_only: bool = False,
) -> str:
    """Build SRA search query for aging-related NGS datasets.
    If nanopore_only is True, restrict to Oxford Nanopore at search time.
    """
    terms = list(AGING_TERMS)
    if extra_terms:
        terms.extend(extra_terms)
    or_part = " OR ".join(f'"{t}"' for t in terms)
    query = f"({or_part}) AND (sequencing OR RNA-Seq OR whole genome OR WGS OR methylation)"
    if nanopore_only:
        query += " AND (nanopore OR \"Oxford Nanopore\")"
    return query


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Search NCBI SRA for aging-related NGS datasets and mark Oxford Nanopore."
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=100,
        help="Maximum number of SRA runs/experiments to fetch (default: 100)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path(__file__).resolve().parent / "aging_ngs_datasets.csv",
        help="Output CSV path for all datasets (default: aging_ngs_datasets.csv)",
    )
    parser.add_argument(
        "--nanopore-output",
        type=Path,
        default=None,
        help="Output CSV path for Oxford Nanopore–only datasets (default: same dir as -o, suffix _nanopore_only.csv)",
    )
    parser.add_argument(
        "--query",
        type=str,
        default=None,
        help="Override search query (default: built-in aging + sequencing terms)",
    )
    args = parser.parse_args()

    # Default nanopore-only path: e.g. aging_ngs_datasets.csv -> aging_ngs_datasets_nanopore_only.csv
    if args.nanopore_output is None:
        stem = args.output.stem
        args.nanopore_output = args.output.parent / f"{stem}_nanopore_only.csv"

    api_key = get_ncbi_api_key()
    if not api_key:
        print("Note: NCBI_API_KEY not set. Consider adding it to .env for higher rate limits.", file=sys.stderr)

    fieldnames = [
        "accession", "title", "description", "link", "number_of_samples",
        "platform", "platform_raw", "is_oxford_nanopore", "organism", "study_accession",
    ]
    header_line = ",".join(fieldnames) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Search 1: all aging-related NGS datasets (no platform filter)
    query_all = args.query or build_search_query(nanopore_only=False)
    print(f"Query (all): {query_all[:100]}...", file=sys.stderr)
    print("Searching SRA (all datasets)...", file=sys.stderr)
    try:
        uid_list_all = search_sra(query_all, max_results=args.max_results, api_key=api_key)
    except requests.RequestException as e:
        print(f"ERROR: SRA search failed: {e}", file=sys.stderr)
        return 1

    # Search 2: same aging terms but restrict to Oxford Nanopore at search time
    if args.query:
        query_nanopore = args.query + ' AND (nanopore OR "Oxford Nanopore")'
    else:
        query_nanopore = build_search_query(nanopore_only=True)
    print(f"Query (Nanopore only): {query_nanopore[:100]}...", file=sys.stderr)
    print("Searching SRA (Oxford Nanopore filter)...", file=sys.stderr)
    time.sleep(0.35)
    try:
        uid_list_nanopore = search_sra(query_nanopore, max_results=args.max_results, api_key=api_key)
    except requests.RequestException as e:
        print(f"ERROR: SRA search (nanopore) failed: {e}", file=sys.stderr)
        return 1

    if not uid_list_all and not uid_list_nanopore:
        print("No results found for either search.", file=sys.stderr)
        with open(args.output, "w", newline="", encoding="utf-8") as f:
            f.write(header_line)
        with open(args.nanopore_output, "w", newline="", encoding="utf-8") as f:
            f.write(header_line)
        return 0

    # Fetch summaries and write CSV for all datasets
    rows_all = []
    if uid_list_all:
        print(f"Found {len(uid_list_all)} IDs (all). Fetching summaries...", file=sys.stderr)
        time.sleep(0.35)
        try:
            rows_all = fetch_sra_summaries(uid_list_all, api_key=api_key)
        except requests.RequestException as e:
            print(f"ERROR: Fetch summaries failed: {e}", file=sys.stderr)
            return 1
        except ET.ParseError as e:
            print(f"ERROR: Failed to parse SRA response: {e}", file=sys.stderr)
            return 1

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows_all)

    # Fetch summaries and write CSV for Nanopore-only (from second search, not filtered in memory)
    rows_nanopore = []
    if uid_list_nanopore:
        print(f"Found {len(uid_list_nanopore)} IDs (Nanopore). Fetching summaries...", file=sys.stderr)
        time.sleep(0.35)
        try:
            rows_nanopore = fetch_sra_summaries(uid_list_nanopore, api_key=api_key)
        except requests.RequestException as e:
            print(f"ERROR: Fetch summaries (nanopore) failed: {e}", file=sys.stderr)
            return 1
        except ET.ParseError as e:
            print(f"ERROR: Failed to parse SRA response (nanopore): {e}", file=sys.stderr)
            return 1

    with open(args.nanopore_output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows_nanopore)

    print(f"Wrote {len(rows_all)} records to {args.output}", file=sys.stderr)
    print(f"Wrote {len(rows_nanopore)} records (Nanopore search) to {args.nanopore_output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
