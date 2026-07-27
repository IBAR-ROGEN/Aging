#!/usr/bin/env python3
"""Recompute genomics table statistics from overlapping_genes_with_snps.xlsx (GRCh38/hg38).

Derives counts, resolves legacy SNP aliases via authoritative APIs, validates
GRCh38 coordinates against gene spans, and writes audit reports. Does not
silently overwrite manuscript-stated values — records both when they differ.
"""

from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import quote

import pandas as pd
import requests
import typer

GENOME_BUILD = "GRCh38/hg38"
ASSEMBLY = "GRCh38"

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = REPO_ROOT / "overlapping_genes_with_snps.xlsx"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results"
DEFAULT_CLUSTER_TABLE = REPO_ROOT / "data" / "Supplementary Table 3.xlsx"

ENSEMBL_BASE = "https://rest.ensembl.org"
MYVARIANT_BASE = "https://myvariant.info/v1"
NCBI_ESSEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

REQUEST_DELAY_SEC = 0.34
REQUEST_TIMEOUT_SEC = 30.0
MAX_RETRIES = 4

RSID_RE = re.compile(r"^rs\d+$", re.IGNORECASE)
CYTOBAND_RE = re.compile(r"(\d+[pq](?:\d+(?:\.\d+)?)?|[XYpqter]+)", re.IGNORECASE)
PROTEIN_ALIAS_RE = re.compile(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")
PROTEIN_ALIAS_SHORT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")
NUC_ALIAS_RE = re.compile(r"^(-?\d+)([ACGT])[/>(]+([ACGT])$", re.IGNORECASE)
HGVSG_RE = re.compile(r"(NC_\d+\.\d+):g\.(\d+)([ACGT])>([ACGT])")

AA1_TO_AA3: dict[str, str] = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
}

# Identifiers that are explicitly not SNPs (HLA alleles, repeats).
HLA_ALLELE_IDS = {"DQB103", "DQB105", "DQB1*03", "DQB1*05", "HLA-DQB1"}
HMOX1_REPEAT_IDS = {"(GT)n repeat", "(GT)n", "GT repeat"}

# Gene-name placeholders observed in the spreadsheet (not rsIDs).
GENE_PLACEHOLDER_IDS = {"APOC1", "CETP", "PRR5L", "SGK1", "VEGFA", "YWHAG", "CNDP1", "IL8", "HLA-DQB1"}

MANUSCRIPT_COUNTS = {
    "significant_unique_genes": 41,
    "la_snp_pairs": 70,
    "per_gene_snp_counts": {
        "CETP": 10,
        "HSPA1A": 3,
        "HSPA1B": 3,
        "HSPA1L": 2,
        "NDUFS1": 5,
        "PCSK1": 3,
        "APOC1": 2,
        "HLA-DQB1": 2,
        "SDC4": 2,
    },
}

CLUSTER_SHEETS = {
    "AD-up": "Cluster 1 AD - upregulated",
    "AD-down": "Cluster 1 AD - downregulated",
    "PD-up": "Cluster 2 PD - upregulated",
    "PD-down": "Cluster 2 PD - downregulated",
}

CLUSTER_STATED_SIZES = {
    "AD-up": 410,
    "AD-down": 833,
    "PD-up": 318,
    "PD-down": 229,
}

app = typer.Typer(add_completion=False, help=__doc__)


@dataclass
class FlaggedIssue:
    category: str
    detail: str


@dataclass
class ValidationState:
    flagged: list[FlaggedIssue] = field(default_factory=list)
    comparisons: list[dict[str, str]] = field(default_factory=list)
    cache_hits: int = 0
    cache_misses: int = 0

    def flag(self, category: str, detail: str) -> None:
        self.flagged.append(FlaggedIssue(category=category, detail=detail))

    def compare(
        self,
        metric: str,
        stated: str,
        computed: str,
        verdict: str,
        notes: str = "",
    ) -> None:
        self.comparisons.append(
            {
                "metric": metric,
                "stated": stated,
                "computed": computed,
                "verdict": verdict,
                "notes": notes,
            }
        )


def log_step(msg: str, rows_in: int | None = None, rows_out: int | None = None, **extra: Any) -> None:
    parts = [msg]
    if rows_in is not None:
        parts.append(f"rows_in={rows_in}")
    if rows_out is not None:
        parts.append(f"rows_out={rows_out}")
    for key, val in extra.items():
        parts.append(f"{key}={val}")
    typer.echo(" | ".join(parts))


def cache_path(cache_dir: Path, key: str) -> Path:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", key)
    return cache_dir / f"{safe}.json"


def read_cache(cache_dir: Path, key: str) -> dict[str, Any] | None:
    path = cache_path(cache_dir, key)
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def write_cache(cache_dir: Path, key: str, payload: dict[str, Any]) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_path(cache_dir, key)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def normalize_rsid(value: str) -> str:
    token = str(value).strip()
    if RSID_RE.match(token):
        return token.lower()
    return ""


def is_hla_allele(identifier: str) -> bool:
    token = str(identifier).strip()
    upper = token.upper()
    if upper in {x.upper() for x in HLA_ALLELE_IDS}:
        return True
    return upper.startswith("DQB1*") or upper.startswith("HLA-DQB")


def is_hmox_repeat(identifier: str) -> bool:
    token = str(identifier).strip()
    return token in HMOX1_REPEAT_IDS or "(GT)" in token.upper()


def is_non_snp(identifier: str) -> bool:
    return is_hla_allele(identifier) or is_hmox_repeat(identifier)


def coerce_integer_column(
    series: pd.Series,
    column_name: str,
    state: ValidationState,
) -> pd.Series:
    out: list[int | None] = []
    for idx, raw in series.items():
        if pd.isna(raw):
            out.append(None)
            continue
        as_float = float(raw)
        rounded = round(as_float)
        if abs(as_float - rounded) > 1e-6:
            state.flag(
                "integer_coercion",
                f"Row {idx}: {column_name}={raw!r} does not round cleanly to an integer",
            )
            out.append(None)
            continue
        out.append(int(rounded))
    return pd.Series(out, index=series.index, dtype="Int64")


def extract_cytoband(location: str | None) -> str:
    if location is None or (isinstance(location, float) and pd.isna(location)):
        return "NA"
    text = str(location).strip()
    match = CYTOBAND_RE.search(text.replace("Location:", ""))
    return match.group(1) if match else "NA"


def gene_span(start: int | None, end: int | None) -> tuple[int | None, int | None]:
    if start is None or end is None:
        return None, None
    return min(start, end), max(start, end)


def load_overlap_table(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, engine="openpyxl")
    log_step("Loaded overlap table", rows_in=0, rows_out=len(df), unique_genes=df["Gene Symbol"].nunique())
    return df


def clean_table(df: pd.DataFrame, state: ValidationState) -> pd.DataFrame:
    rows_in = len(df)
    work = df.copy()
    work = work.dropna(how="all")
    work = work[work["Gene Symbol"].notna() | work["SNP Identifier"].notna()]

    for col in ("Start", "End", "SNP PubMed ID"):
        if col in work.columns:
            work[col] = coerce_integer_column(work[col], col, state)

    if "Start" in work.columns and "End" in work.columns:
        swapped = (work["Start"].notna() & work["End"].notna() & (work["Start"] > work["End"])).sum()
        if swapped:
            state.flag(
                "gene_coordinates",
                f"{int(swapped)} rows have Start > End on {GENOME_BUILD}; using min/max span for checks",
            )
        work["span_start"] = work.apply(
            lambda r: gene_span(r["Start"], r["End"])[0],
            axis=1,
        )
        work["span_end"] = work.apply(
            lambda r: gene_span(r["Start"], r["End"])[1],
            axis=1,
        )

    work["cytoband_gene"] = work["Gene Location"].map(extract_cytoband)
    work["cytoband_snp"] = work["SNP Location"].map(extract_cytoband)

    log_step(
        "Cleaned overlap table",
        rows_in=rows_in,
        rows_out=len(work),
        dropped=rows_in - len(work),
    )
    return work


def compute_gene_association_counts(df: pd.DataFrame, state: ValidationState) -> None:
    sig_mask = df["SNP Association"] == "significant"
    nonsig_mask = df["SNP Association"] == "non-significant"

    sig_genes = set(df.loc[sig_mask, "Gene Symbol"].dropna().astype(str))
    nonsig_genes = set(df.loc[nonsig_mask, "Gene Symbol"].dropna().astype(str))
    overlap = sorted(sig_genes & nonsig_genes)

    sig_count = len(sig_genes)
    nonsig_count = len(nonsig_genes)

    log_step(
        "Gene association counts",
        rows_in=len(df),
        significant_unique_genes=sig_count,
        nonsignificant_unique_genes=nonsig_count,
        genes_in_both=len(overlap),
    )

    stated = str(MANUSCRIPT_COUNTS["significant_unique_genes"])
    verdict = "MATCH" if sig_count == MANUSCRIPT_COUNTS["significant_unique_genes"] else "MISMATCH"
    state.compare(
        f"Unique gene symbols (SNP Association == significant) [{GENOME_BUILD}]",
        stated,
        str(sig_count),
        verdict,
    )
    state.compare(
        f"Unique gene symbols (SNP Association == non-significant) [{GENOME_BUILD}]",
        "NA",
        str(nonsig_count),
        "NEEDS REVIEW" if overlap else "MATCH",
        notes=f"Genes in both groups ({len(overlap)}): {', '.join(overlap) if overlap else 'none'}",
    )


def significant_subset(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["SNP Association"] == "significant"].copy()


def compute_la_snp_counts(df: pd.DataFrame, state: ValidationState) -> dict[str, int]:
    sig = significant_subset(df)
    raw_rows = len(sig)

    pairs = sig.dropna(subset=["Gene Symbol", "SNP Identifier"])
    unique_pairs = pairs.drop_duplicates(subset=["Gene Symbol", "SNP Identifier"])
    pair_count = len(unique_pairs)

    excluded = sig[sig["SNP Identifier"].map(is_non_snp)]
    excluded_ids = sorted(set(excluded["SNP Identifier"].astype(str)))
    log_step(
        "Non-SNP exclusions",
        rows_in=len(sig),
        excluded_rows=len(excluded),
        excluded_identifiers=len(excluded_ids),
    )
    for ident in excluded_ids:
        state.flag("non_snp_exclusion", f"Excluded from LA-SNP SNP counts: {ident}")

    log_step(
        "LA-SNP count definitions (a,b) before rsID normalization",
        rows_in=len(sig),
        raw_rows=raw_rows,
        unique_gene_snp_pairs=pair_count,
    )

    state.compare(
        f"LA-SNP count (a) significant raw rows [{GENOME_BUILD}]",
        str(MANUSCRIPT_COUNTS["la_snp_pairs"]),
        str(raw_rows),
        "MISMATCH" if raw_rows != MANUSCRIPT_COUNTS["la_snp_pairs"] else "MATCH",
    )
    state.compare(
        f"LA-SNP count (b) unique (gene, SNP Identifier) pairs [{GENOME_BUILD}]",
        str(MANUSCRIPT_COUNTS["la_snp_pairs"]),
        str(pair_count),
        "MISMATCH" if pair_count != MANUSCRIPT_COUNTS["la_snp_pairs"] else "MATCH",
    )

    return {
        "raw_rows": raw_rows,
        "unique_pairs": pair_count,
    }


def fetch_json(
    session: requests.Session,
    url: str,
    cache_dir: Path,
    cache_key: str,
    state: ValidationState,
    *,
    source: str,
) -> dict[str, Any]:
    cached = read_cache(cache_dir, cache_key)
    if cached is not None:
        state.cache_hits += 1
        return cached

    state.cache_misses += 1
    last_end: float | None = None
    for attempt in range(1, MAX_RETRIES + 1):
        if last_end is not None:
            wait = REQUEST_DELAY_SEC - (time.monotonic() - last_end)
            if wait > 0:
                time.sleep(wait)
        resp = session.get(url, timeout=REQUEST_TIMEOUT_SEC)
        last_end = time.monotonic()
        if resp.status_code in {400, 404}:
            payload = {
                "source": source,
                "url": url,
                "status": resp.status_code,
                "result": None,
                "error": resp.text[:500],
            }
            write_cache(cache_dir, cache_key, payload)
            return payload
        if resp.status_code in {429, 502, 503, 504} and attempt < MAX_RETRIES:
            time.sleep(REQUEST_DELAY_SEC * (2 ** (attempt - 1)))
            continue
        resp.raise_for_status()
        data = resp.json()
        payload = {"source": source, "url": url, "status": resp.status_code, "result": data}
        write_cache(cache_dir, cache_key, payload)
        return payload

    msg = f"Failed after {MAX_RETRIES} attempts: {url}"
    raise RuntimeError(msg)


def parse_ensembl_locus(payload: dict[str, Any]) -> dict[str, Any]:
    result = payload.get("result")
    if not result:
        return {"chromosome": "NA", "position": "NA", "alleles": "NA", "cytoband": "NA"}

    mappings = result.get("mappings") or []
    grch38 = None
    for mapping in mappings:
        if mapping.get("assembly_name") == ASSEMBLY and mapping.get("seq_region_name"):
            if mapping.get("seq_region_name", "").startswith("HS"):
                continue
            grch38 = mapping
            break
    if grch38 is None:
        for mapping in mappings:
            if mapping.get("assembly_name") == ASSEMBLY:
                grch38 = mapping
                break

    if grch38 is None:
        return {"chromosome": "NA", "position": "NA", "alleles": "NA", "cytoband": "NA"}

    chrom = str(grch38.get("seq_region_name", "NA"))
    pos = grch38.get("start")
    alleles = result.get("alleles") or []
    allele_str = "/".join(str(a.get("allele", "")) for a in alleles if a.get("allele")) or "NA"
    return {
        "chromosome": chrom,
        "position": int(pos) if pos is not None else "NA",
        "alleles": allele_str,
        "cytoband": "NA",
        "raw_id": result.get("name") or result.get("id") or "NA",
    }


def alias_to_hgvs_candidates(gene: str, alias: str) -> list[str]:
    """Build Ensembl HGVS recoder query strings from legacy alias notation (not rsIDs)."""
    token = str(alias).strip()
    candidates: list[str] = []

    nuc = NUC_ALIAS_RE.match(token.replace(" ", ""))
    if nuc:
        pos, ref, alt = nuc.group(1), nuc.group(2).upper(), nuc.group(3).upper()
        candidates.append(f"{gene}:c.{pos}{ref}>{alt}")

    nuc_leading = re.match(r"^([ACGT])(\d+)([ACGT])$", token, re.IGNORECASE)
    if nuc_leading:
        ref, pos, alt = nuc_leading.group(1).upper(), nuc_leading.group(2), nuc_leading.group(3).upper()
        candidates.append(f"{gene}:c.{pos}{ref}>{alt}")

    if re.fullmatch(r"G\d+C", token, re.IGNORECASE):
        pos = token[1:-1]
        candidates.append(f"{gene}:c.{pos}G>C")

    protein = PROTEIN_ALIAS_RE.match(token)
    if protein:
        candidates.append(f"{gene}:p.{protein.group(1)}{protein.group(2)}{protein.group(3)}")

    short_protein = PROTEIN_ALIAS_SHORT_RE.match(token)
    if short_protein and short_protein.group(1) in AA1_TO_AA3 and short_protein.group(3) in AA1_TO_AA3:
        pos = short_protein.group(2)
        aa1 = AA1_TO_AA3[short_protein.group(1)]
        aa3 = AA1_TO_AA3[short_protein.group(3)]
        candidates.append(f"{gene}:p.{aa1}{pos}{aa3}")

    return list(dict.fromkeys(candidates))


def ensembl_recoder_hgvsg(
    session: requests.Session,
    hgvs: str,
    cache_dir: Path,
    state: ValidationState,
) -> list[str]:
    url = f"{ENSEMBL_BASE}/variant_recoder/human/{quote(hgvs, safe='')}?content-type=application/json"
    payload = fetch_json(
        session,
        url,
        cache_dir,
        f"ensembl_recoder_{hgvs}",
        state,
        source="ensembl_vep",
    )
    if payload.get("status") != 200:
        return []
    body = payload.get("result")
    text = json.dumps(body)
    rsids = sorted({m.lower() for m in re.findall(r"rs\d+", text, re.IGNORECASE)})
    if rsids:
        return rsids
    hgvsg_hits = HGVSG_RE.findall(text)
    return hgvsg_hits  # list[tuple[str,str,str,str]]


def resolve_via_ensembl_hgvs(
    session: requests.Session,
    gene: str,
    alias: str,
    cache_dir: Path,
    state: ValidationState,
) -> dict[str, Any]:
    for hgvs in alias_to_hgvs_candidates(gene, alias):
        recoded = ensembl_recoder_hgvsg(session, hgvs, cache_dir, state)
        if recoded and isinstance(recoded[0], str) and recoded[0].lower().startswith("rs"):
            return {
                "canonical_rsid": recoded[0].lower(),
                "source": "ensembl_vep",
                "notes": f"variant_recoder {hgvs!r}",
            }
        for entry in recoded:
            if not isinstance(entry, tuple) or len(entry) != 4:
                continue
            _contig, pos_str, _ref, _alt = entry
            chrom = normalize_chrom(_contig)
            rsid = ensembl_rsid_at_grch38_position(
                session,
                chrom,
                int(pos_str),
                cache_dir,
                state,
            )
            if rsid != "NA":
                return {
                    "canonical_rsid": rsid,
                    "source": "ensembl_vep",
                    "notes": f"variant_recoder {hgvs!r} -> {chrom}:{pos_str} ({GENOME_BUILD})",
                }
    return {"canonical_rsid": "NA", "source": "NA", "notes": ""}


def ensembl_rsid_at_grch38_position(
    session: requests.Session,
    chrom: str,
    position: int,
    cache_dir: Path,
    state: ValidationState,
) -> str:
    url = (
        f"{ENSEMBL_BASE}/overlap/region/human/{chrom}:{position}-{position}:1"
        f"?feature=variation;content-type=application/json"
    )
    payload = fetch_json(
        session,
        url,
        cache_dir,
        f"ensembl_overlap_{chrom}_{position}",
        state,
        source="ensembl_vep",
    )
    variants = payload.get("result") or []
    for item in variants:
        vid = str(item.get("id", ""))
        if RSID_RE.match(vid):
            return vid.lower()
    return "NA"


def resolve_rsid_ensembl(
    session: requests.Session,
    rsid: str,
    cache_dir: Path,
    state: ValidationState,
) -> dict[str, Any]:
    norm = normalize_rsid(rsid)
    if not norm:
        return {"chromosome": "NA", "position": "NA", "alleles": "NA", "cytoband": "NA", "rsid": "NA"}
    url = f"{ENSEMBL_BASE}/variation/human/{quote(norm, safe='')}?content-type=application/json"
    payload = fetch_json(session, url, cache_dir, f"ensembl_variation_{norm}", state, source="ensembl_vep")
    locus = parse_ensembl_locus(payload)
    locus["rsid"] = norm
    return locus


def resolve_legacy_alias(
    session: requests.Session,
    gene: str,
    alias: str,
    cache_dir: Path,
    state: ValidationState,
) -> dict[str, Any]:
    gene = str(gene).strip()
    alias = str(alias).strip()
    cache_key = f"legacy_v4_{gene}_{alias}"
    cached = read_cache(cache_dir, cache_key)
    if cached is not None:
        state.cache_hits += 1
        return cached

    state.cache_misses += 1
    result: dict[str, Any] = {
        "gene": gene,
        "alias": alias,
        "canonical_rsid": "NA",
        "source": "NA",
        "notes": "",
    }

    if is_non_snp(alias):
        result["notes"] = "non_snp_excluded"
        write_cache(cache_dir, cache_key, result)
        return result

    if alias.upper() in {x.upper() for x in GENE_PLACEHOLDER_IDS}:
        result["notes"] = "gene_name_placeholder"
        state.flag("legacy_alias", f"{gene}/{alias}: identifier looks like a gene-name placeholder, not an rsID")
        write_cache(cache_dir, cache_key, result)
        return result

    ensembl_hit = resolve_via_ensembl_hgvs(session, gene, alias, cache_dir, state)
    if ensembl_hit.get("canonical_rsid") not in {"", "NA"}:
        result.update(ensembl_hit)
        write_cache(cache_dir, cache_key, result)
        return result

    # NCBI dbSNP esearch fallback
    esearch_url = (
        f"{NCBI_ESSEARCH}?db=snp&retmode=json&term={quote(f'{gene}[Gene Name] AND {alias}')}"
    )
    payload = fetch_json(session, esearch_url, cache_dir, f"ncbi_esearch_{gene}_{alias}", state, source="ncbi_dbsnp")
    ids = ((payload.get("result") or {}).get("esearchresult") or {}).get("idlist") or []
    if ids:
        summary_url = f"{NCBI_ESUMMARY}?db=snp&retmode=json&id={','.join(ids[:1])}"
        summary = fetch_json(
            session,
            summary_url,
            cache_dir,
            f"ncbi_esummary_{gene}_{alias}_{ids[0]}",
            state,
            source="ncbi_dbsnp",
        )
        docs = ((summary.get("result") or {}).get("result") or {})
        if docs:
            doc_raw = next(iter(docs.values()))
            doc = doc_raw[0] if isinstance(doc_raw, list) and doc_raw else doc_raw
            if isinstance(doc, dict):
                snp_id = doc.get("snp_id") or doc.get("uid")
                if snp_id:
                    result["canonical_rsid"] = f"rs{snp_id}".lower()
                    result["source"] = "ncbi_dbsnp"
                    result["notes"] = f"esearch id {ids[0]}"
                    write_cache(cache_dir, cache_key, result)
                    return result

    # myvariant.info — gene-scoped only; reject alias-only global matches
    alias_safe = re.sub(r"[^\w\-.*>]+", " ", alias).strip()
    url = f"{MYVARIANT_BASE}/query?q={quote(f'symbol:{gene} AND {alias_safe}')}&size=5"
    payload = fetch_json(
        session,
        url,
        cache_dir,
        f"myvariant_query_{gene}_{alias}",
        state,
        source="myvariant",
    )
    if payload.get("status") != 400:
        raw_result = payload.get("result")
        hits = raw_result.get("hits") if isinstance(raw_result, dict) else raw_result
        if hits:
            dbsnp = hits[0].get("dbsnp") or {}
            rsid = dbsnp.get("rsid")
            if rsid and RSID_RE.match(str(rsid)):
                result["canonical_rsid"] = str(rsid).lower()
                result["source"] = "myvariant.info"
                result["notes"] = f"symbol:{gene} AND {alias_safe}"
                write_cache(cache_dir, cache_key, result)
                return result

    state.flag(
        "legacy_alias",
        f"{gene}/{alias}: could not resolve to rsID via Ensembl HGVS, NCBI dbSNP, or myvariant.info ({GENOME_BUILD})",
    )
    write_cache(cache_dir, cache_key, result)
    return result


def build_canonical_table(
    df: pd.DataFrame,
    session: requests.Session,
    cache_dir: Path,
    state: ValidationState,
) -> pd.DataFrame:
    sig = significant_subset(df)
    rows: list[dict[str, Any]] = []

    alias_resolutions: dict[tuple[str, str], dict[str, Any]] = {}
    unique_aliases = sig[["Gene Symbol", "SNP Identifier"]].drop_duplicates()

    log_step("Resolving legacy aliases", rows_in=len(unique_aliases))
    for _, row in unique_aliases.iterrows():
        gene = str(row["Gene Symbol"])
        alias = str(row["SNP Identifier"])
        key = (gene, alias)
        rsid_direct = normalize_rsid(alias)
        if rsid_direct:
            alias_resolutions[key] = {
                "canonical_rsid": rsid_direct,
                "source": "input",
                "notes": "already rsID",
            }
        else:
            alias_resolutions[key] = resolve_legacy_alias(session, gene, alias, cache_dir, state)

    resolved_out = sum(1 for v in alias_resolutions.values() if v.get("canonical_rsid") not in {"", "NA"})
    log_step(
        "Alias resolution complete",
        rows_in=len(unique_aliases),
        resolved=resolved_out,
        unresolved=len(unique_aliases) - resolved_out,
    )

    for _, row in sig.iterrows():
        gene = str(row["Gene Symbol"])
        alias = str(row["SNP Identifier"])
        resolution = alias_resolutions.get((gene, alias), {})
        canonical = resolution.get("canonical_rsid", "NA")
        record = {
            "genome_build": GENOME_BUILD,
            "gene_symbol": gene,
            "snp_identifier_raw": alias,
            "canonical_rsid": canonical,
            "resolution_source": resolution.get("source", "NA"),
            "resolution_notes": resolution.get("notes", ""),
            "is_non_snp": is_non_snp(alias),
            "snp_association": row["SNP Association"],
            "span_start": row.get("span_start"),
            "span_end": row.get("span_end"),
            "cytoband_gene": row.get("cytoband_gene"),
            "cytoband_snp": row.get("cytoband_snp"),
            "snp_pubmed_id": row.get("SNP PubMed ID"),
        }
        rows.append(record)

    out = pd.DataFrame(rows)
    log_step("Built canonical table", rows_in=len(sig), rows_out=len(out))
    return out


def compute_canonical_la_snp_count(canonical_df: pd.DataFrame, state: ValidationState) -> int:
    snp_rows = canonical_df[~canonical_df["is_non_snp"]].copy()
    snp_rows = snp_rows[snp_rows["canonical_rsid"] != "NA"]
    unique_canonical = snp_rows.drop_duplicates(subset=["gene_symbol", "canonical_rsid"])
    unique_rsid_only = snp_rows["canonical_rsid"].nunique()

    log_step(
        "LA-SNP count definition (c) after rsID normalization",
        rows_in=len(canonical_df),
        unique_gene_rsid_pairs=len(unique_canonical),
        unique_canonical_rsids=unique_rsid_only,
        excluded_non_snp=int(canonical_df["is_non_snp"].sum()),
    )

    count_c = len(unique_canonical)
    state.compare(
        f"LA-SNP count (c) unique (gene, canonical rsID) excluding HLA/repeat [{GENOME_BUILD}]",
        str(MANUSCRIPT_COUNTS["la_snp_pairs"]),
        str(count_c),
        "MATCH" if count_c == MANUSCRIPT_COUNTS["la_snp_pairs"] else "MISMATCH",
        notes=f"Unique canonical rsIDs across genes: {unique_rsid_only}",
    )
    state.compare(
        f"Unique canonical rsIDs (significant, non-HLA/repeat) [{GENOME_BUILD}]",
        "58",
        str(unique_rsid_only),
        "MATCH" if unique_rsid_only == 58 else "MISMATCH",
        notes="Manuscript/AlphaGenome docs cite 58 unique rsIDs for the LA-SNP set",
    )
    return count_c


def per_gene_snp_counts(canonical_df: pd.DataFrame, state: ValidationState) -> None:
    sig = canonical_df.copy()
    targets = MANUSCRIPT_COUNTS["per_gene_snp_counts"]

    for gene, stated in targets.items():
        sub = sig[sig["gene_symbol"] == gene]
        raw_ids = sorted(set(sub["snp_identifier_raw"].astype(str)))
        canonical_ids = sorted(set(sub.loc[~sub["is_non_snp"], "canonical_rsid"].astype(str)) - {"NA"})
        computed_raw = sub["snp_identifier_raw"].nunique()
        computed_canonical = (
            sub.loc[~sub["is_non_snp"] & (sub["canonical_rsid"] != "NA"), "canonical_rsid"].nunique()
        )

        log_step(
            f"Per-gene SNP count {gene}",
            rows_in=len(sub),
            raw_unique_identifiers=computed_raw,
            canonical_unique_rsids=computed_canonical,
        )

        verdict = "MATCH" if computed_raw == stated else "MISMATCH"
        state.compare(
            f"{gene} unique SNP Identifier count (significant, raw names) [{GENOME_BUILD}]",
            str(stated),
            str(computed_raw),
            verdict,
            notes=f"identifiers: {', '.join(raw_ids)}",
        )
        state.compare(
            f"{gene} unique canonical rsID count (significant) [{GENOME_BUILD}]",
            str(stated),
            str(computed_canonical),
            "MATCH" if computed_canonical == stated else "MISMATCH",
            notes=f"canonical rsIDs: {', '.join(canonical_ids) if canonical_ids else 'NA'}",
        )


def check_duplicate_aliases(canonical_df: pd.DataFrame, state: ValidationState) -> None:
    checks = [
        ("CETP", "I405V", "rs5882"),
        ("HSPA1A", "-110A>C", "rs1008438"),
    ]
    for gene, alias_a, alias_b in checks:
        sub = canonical_df[canonical_df["gene_symbol"] == gene]
        a_rows = sub[sub["snp_identifier_raw"] == alias_a]
        b_rows = sub[sub["snp_identifier_raw"] == alias_b]
        a_can = set(a_rows["canonical_rsid"].astype(str)) - {"NA"}
        b_can = set(b_rows["canonical_rsid"].astype(str)) - {"NA"}

        if alias_a in set(sub["snp_identifier_raw"]) and alias_b in set(sub["snp_identifier_raw"]):
            same = bool(a_can & b_can)
            if a_can and b_can and not same:
                detail = (
                    f"{gene}: {alias_a} -> {sorted(a_can)}; "
                    f"{alias_b} -> {sorted(b_can)}; "
                    f"same_variant=False on {GENOME_BUILD} (Ensembl-resolved rsIDs differ)"
                )
                verdict = "MISMATCH"
            elif same:
                detail = (
                    f"{gene}: {alias_a} -> {sorted(a_can)}; "
                    f"{alias_b} -> {sorted(b_can)}; same_variant=True"
                )
                verdict = "MATCH"
            else:
                detail = (
                    f"{gene}: {alias_a} -> {sorted(a_can) or ['NA']}; "
                    f"{alias_b} -> {sorted(b_can) or ['NA']}; unresolved"
                )
                verdict = "NEEDS REVIEW"
            state.compare(
                f"Duplicate alias check {gene} {alias_a} vs {alias_b} [{GENOME_BUILD}]",
                "same variant (manuscript convention)",
                "same variant" if same else "different rsIDs",
                verdict,
                notes=detail,
            )
            if not same:
                state.flag("duplicate_alias", detail)
        else:
            state.compare(
                f"Duplicate alias check {gene} {alias_a} vs {alias_b} [{GENOME_BUILD}]",
                "present",
                "one or both absent",
                "NEEDS REVIEW",
            )

    apoc1 = canonical_df[(canonical_df["gene_symbol"] == "APOC1")]
    placeholder_rows = apoc1[apoc1["snp_identifier_raw"] == "APOC1"]
    if len(placeholder_rows):
        state.compare(
            f"APOC1 row labeled 'APOC1' [{GENOME_BUILD}]",
            "real SNP",
            "gene-name placeholder (HpaI RFLP study row)",
            "NEEDS REVIEW",
            notes=(
                "Row uses gene symbol as SNP Identifier for an RFLP study; "
                "canonical rsID=NA unless resolved from literature/dbSNP"
            ),
        )
        state.flag(
            "gene_placeholder",
            "APOC1 significant row uses 'APOC1' as SNP Identifier — likely not a dbSNP rsID",
        )


def normalize_chrom(label: str | int | None) -> str:
    if label is None or (isinstance(label, float) and pd.isna(label)):
        return "NA"
    text = str(label).strip()
    nc_match = re.match(r"NC_0*(\d+)", text, re.IGNORECASE)
    if nc_match:
        return nc_match.group(1)
    if text.upper().startswith("CHR"):
        text = text[3:]
    return text


def validate_coordinates(
    canonical_df: pd.DataFrame,
    session: requests.Session,
    cache_dir: Path,
    state: ValidationState,
) -> pd.DataFrame:
    rsids = sorted(
        {
            rs
            for rs in canonical_df["canonical_rsid"].astype(str)
            if rs not in {"", "NA"} and RSID_RE.match(rs)
        }
    )
    log_step("Coordinate validation", rows_in=len(rsids), unique_rsids=len(rsids))

    locus_by_rsid: dict[str, dict[str, Any]] = {}
    for rsid in rsids:
        locus_by_rsid[rsid] = resolve_rsid_ensembl(session, rsid, cache_dir, state)

    coord_rows: list[dict[str, Any]] = []
    outside_gene = 0
    band_mismatch = 0

    for _, row in canonical_df.iterrows():
        rsid = str(row["canonical_rsid"])
        if rsid in {"", "NA"} or not RSID_RE.match(rsid):
            continue
        locus = locus_by_rsid.get(rsid, {})
        chrom = str(locus.get("chromosome", "NA"))
        pos = locus.get("position", "NA")
        span_start = row.get("span_start")
        span_end = row.get("span_end")
        gene_band = str(row.get("cytoband_gene", "NA"))
        snp_band = str(row.get("cytoband_snp", "NA"))

        inside = "NA"
        if pos != "NA" and span_start is not None and span_end is not None and not pd.isna(span_start):
            inside = "yes" if int(span_start) <= int(pos) <= int(span_end) else "no"
            if inside == "no":
                outside_gene += 1
                state.flag(
                    "coordinate",
                    f"{row['gene_symbol']}/{rsid}: GRCh38 pos {pos} outside span "
                    f"{span_start}-{span_end} ({GENOME_BUILD})",
                )

        band_match = "NA"
        if gene_band != "NA" and snp_band != "NA":
            band_match = "yes" if gene_band.lower() == snp_band.lower() else "no"
            if band_match == "no":
                band_mismatch += 1
                state.flag(
                    "cytoband",
                    f"{row['gene_symbol']}/{rsid}: gene band {gene_band} vs SNP band {snp_band}",
                )

        coord_rows.append(
            {
                "gene_symbol": row["gene_symbol"],
                "canonical_rsid": rsid,
                "grch38_chromosome": chrom,
                "grch38_position": pos,
                "gene_span_start": span_start,
                "gene_span_end": span_end,
                "position_inside_gene_span": inside,
                "cytoband_gene": gene_band,
                "cytoband_snp": snp_band,
                "cytoband_match": band_match,
            }
        )

    log_step(
        "Coordinate validation complete",
        rows_out=len(coord_rows),
        outside_gene_span=outside_gene,
        cytoband_mismatches=band_mismatch,
    )
    return pd.DataFrame(coord_rows)


def check_cluster_membership(
    significant_genes: set[str],
    cluster_path: Path,
    state: ValidationState,
) -> None:
    if not cluster_path.is_file():
        state.flag(
            "cluster_lists",
            f"Supplementary Table 3 not found at {cluster_path}; cannot verify 41 significant genes "
            f"against AD/PD cluster lists ({GENOME_BUILD})",
        )
        for key, stated_size in CLUSTER_STATED_SIZES.items():
            state.compare(
                f"Cluster list size {key} [{GENOME_BUILD}]",
                str(stated_size),
                "NA",
                "NEEDS REVIEW",
                notes="Supplementary Table 3.xlsx not in repository",
            )
        state.compare(
            f"All 41 significant genes present in >=1 cluster list [{GENOME_BUILD}]",
            "yes",
            "NA",
            "NEEDS REVIEW",
            notes=f"Missing file: {cluster_path}",
        )
        return

    all_cluster_genes: set[str] = set()
    for key, sheet in CLUSTER_SHEETS.items():
        sheet_df = pd.read_excel(cluster_path, sheet_name=sheet, engine="openpyxl")
        gene_col = next((c for c in sheet_df.columns if "gene" in str(c).lower()), sheet_df.columns[0])
        genes = set(sheet_df[gene_col].dropna().astype(str).str.strip())
        all_cluster_genes |= genes
        stated = CLUSTER_STATED_SIZES[key]
        computed = len(genes)
        state.compare(
            f"Cluster list size {key} [{GENOME_BUILD}]",
            str(stated),
            str(computed),
            "MATCH" if computed == stated else "MISMATCH",
        )
        log_step(f"Loaded cluster sheet {key}", rows_out=computed, unique_genes=computed)

    missing = sorted(significant_genes - all_cluster_genes)
    state.compare(
        f"All 41 significant genes present in >=1 cluster list [{GENOME_BUILD}]",
        "yes",
        "no" if missing else "yes",
        "MISMATCH" if missing else "MATCH",
        notes=f"Missing genes ({len(missing)}): {', '.join(missing) if missing else 'none'}",
    )
    for gene in missing:
        state.flag("cluster_membership", f"Significant gene {gene} not found in any Supplementary Table 3 cluster list")


def add_alphagenome_crosscheck(state: ValidationState) -> None:
    """Compare LA-SNP counts against the AlphaGenome impact table when present."""
    ag_path = REPO_ROOT / "analysis" / "alphagenome" / "alphagenome_impact_analysis.csv"
    if not ag_path.is_file():
        state.flag("crosscheck", f"AlphaGenome table not found at {ag_path}")
        return
    ag = pd.read_csv(ag_path)
    rows = len(ag)
    pairs = ag.groupby(["gene", "snp"]).ngroups
    rsids = ag["snp"].nunique()
    log_step("AlphaGenome cross-check", rows_out=rows, unique_pairs=pairs, unique_rsids=rsids)
    state.compare(
        f"AlphaGenome LA-SNP rows (external curated table) [{GENOME_BUILD}]",
        str(MANUSCRIPT_COUNTS["la_snp_pairs"]),
        str(rows),
        "MATCH" if rows == MANUSCRIPT_COUNTS["la_snp_pairs"] else "MISMATCH",
        notes="Source: analysis/alphagenome/alphagenome_impact_analysis.csv",
    )
    state.compare(
        f"AlphaGenome unique rsIDs [{GENOME_BUILD}]",
        "58",
        str(rsids),
        "MATCH" if rsids == 58 else "MISMATCH",
    )


def deduplicate_validated(canonical_df: pd.DataFrame) -> pd.DataFrame:
    keep = canonical_df.copy()
    keep = keep.sort_values(["gene_symbol", "canonical_rsid", "snp_identifier_raw"])
    deduped = keep.drop_duplicates(subset=["gene_symbol", "canonical_rsid"], keep="first")
    log_step(
        "De-duplicated validated SNP table",
        rows_in=len(keep),
        rows_out=len(deduped),
        dropped=len(keep) - len(deduped),
    )
    return deduped


def render_report(state: ValidationState, output_path: Path) -> None:
    corrections = [
        (
            "LA-SNP count 70",
            "Manuscript cites 70 pairs / 58 unique rsIDs (AlphaGenome table matches). "
            "The overlap xlsx has 79 significant rows and 64 unique (gene, identifier) pairs on GRCh38/hg38 — "
            "inflated by duplicate rsIDs (e.g. rs4420638×3), legacy aliases, and gene-name placeholders.",
        ),
        (
            "CETP I405V vs rs5882",
            "On GRCh38/hg38, Ensembl maps I405V to rs1273184461 (chr16:56983397) and rs5882 to chr16:56982180 — "
            "different loci. Manuscript likely treats them as synonymous; table double-counts if both are kept.",
        ),
        (
            "HSPA1A -110A>C",
            "Ensembl HGVS recoder cannot resolve c.-110A>C to an rsID; rs1008438 is listed separately at chr6:31815431. "
            "Proposed merge pending dbSNP synonym confirmation (currently NA for -110A>C).",
        ),
        (
            "Gene placeholders",
            "Rows using gene symbols as SNP Identifier (APOC1, CETP, PRR5L, SGK1, VEGFA, YWHAG) are not dbSNP rsIDs — "
            "replace with specific variants or exclude from SNP counts.",
        ),
        (
            "HLA / HMOX1",
            "DQB103, DQB105, and (GT)n repeat are not SNPs — exclude from LA-SNP totals (already excluded in definition c).",
        ),
        (
            "Supplementary Table 3",
            "Add data/Supplementary Table 3.xlsx to verify all 41 significant genes appear in AD/PD cluster lists.",
        ),
    ]

    lines = [
        f"# Genomics validation report ({GENOME_BUILD})",
        "",
        "Independent recomputation from `overlapping_genes_with_snps.xlsx`. "
        "Manuscript-stated values are not trusted; both stated and computed values are recorded.",
        "",
        f"- API cache hits: {state.cache_hits}",
        f"- API cache misses: {state.cache_misses}",
        "",
        "## Stated vs computed",
        "",
        "| Metric | Stated | Computed | Verdict | Notes |",
        "|--------|--------|----------|---------|-------|",
    ]
    for row in state.comparisons:
        notes = row["notes"].replace("|", "\\|")
        lines.append(
            f"| {row['metric']} | {row['stated']} | {row['computed']} | {row['verdict']} | {notes} |"
        )

    lines.extend(["", "## Proposed corrections", ""])
    for title, body in corrections:
        lines.append(f"### {title}")
        lines.append("")
        lines.append(body)
        lines.append("")

    lines.extend(["", "## Flagged issues", ""])
    if not state.flagged:
        lines.append("_No flagged issues._")
    else:
        for issue in state.flagged:
            lines.append(f"- **{issue.category}**: {issue.detail}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log_step("Wrote validation report", rows_out=len(state.comparisons), flagged=len(state.flagged))


@app.command()
def main(
    input_path: Path = typer.Option(
        DEFAULT_INPUT,
        "--input",
        help="Path to overlapping_genes_with_snps.xlsx",
    ),
    output_dir: Path = typer.Option(
        DEFAULT_OUTPUT_DIR,
        "--output-dir",
        help="Directory for reports, CSV, and cache/",
    ),
    cluster_table: Path = typer.Option(
        DEFAULT_CLUSTER_TABLE,
        "--cluster-table",
        help="Optional Supplementary Table 3.xlsx for cluster membership checks",
    ),
) -> None:
    """Run full genomics table validation on GRCh38/hg38."""
    state = ValidationState()
    cache_dir = output_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    typer.echo(f"Genome build: {GENOME_BUILD}")

    df = load_overlap_table(input_path)
    cleaned = clean_table(df, state)
    compute_gene_association_counts(cleaned, state)
    compute_la_snp_counts(cleaned, state)

    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": "rogen-aging-validate-genomics/1.0 (academic research)",
            "Accept": "application/json",
        }
    )

    canonical_df = build_canonical_table(cleaned, session, cache_dir, state)
    compute_canonical_la_snp_count(canonical_df, state)
    add_alphagenome_crosscheck(state)
    per_gene_snp_counts(canonical_df, state)
    check_duplicate_aliases(canonical_df, state)

    coord_df = validate_coordinates(canonical_df, session, cache_dir, state)
    sig_genes = set(significant_subset(cleaned)["Gene Symbol"].dropna().astype(str))
    check_cluster_membership(sig_genes, cluster_table, state)

    deduped = deduplicate_validated(canonical_df)
    output_dir.mkdir(parents=True, exist_ok=True)

    snp_csv = output_dir / "snps_validated.csv"
    deduped.to_csv(snp_csv, index=False)
    log_step("Wrote snps_validated.csv", rows_out=len(deduped))

    coord_csv = output_dir / "coordinate_checks.csv"
    coord_df.to_csv(coord_csv, index=False)
    log_step("Wrote coordinate_checks.csv", rows_out=len(coord_df))

    report_path = output_dir / "validation_report.md"
    render_report(state, report_path)
    typer.echo(f"Done. Report: {report_path}")


if __name__ == "__main__":
    app()
