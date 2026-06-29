#!/usr/bin/env python3
"""Annotate validated longevity SNPs: VEP consequences, AlphaMissense, GWAS, coding vs regulatory split."""

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
REPO_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_INPUT = REPO_ROOT / "results" / "snps_validated.csv"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results"
DEFAULT_CACHE_DIR = DEFAULT_OUTPUT_DIR / "cache" / "variant_annotation"
LEGACY_VEP_CACHE = REPO_ROOT / "analysis" / "vep_cache"

ENSEMBL_BASE = "https://rest.ensembl.org"
MYVARIANT_BASE = "https://myvariant.info/v1"
GWAS_BASE = "https://www.ebi.ac.uk/gwas/rest/api"

REQUEST_DELAY_SEC = 0.34
REQUEST_TIMEOUT_SEC = 30.0
MAX_RETRIES = 4

RSID_RE = re.compile(r"^rs\d+$", re.IGNORECASE)

IMPACT_RANK: dict[str, int] = {
    "HIGH": 0,
    "MODERATE": 1,
    "LOW": 2,
    "MODIFIER": 3,
}

MISSENSE_TERMS = {"missense_variant"}
SYNONYMOUS_TERMS = {"synonymous_variant"}
UTR_TERMS = {
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "5_prime_UTR_truncation",
    "3_prime_UTR_truncation",
}
INTRONIC_TERMS = {"intron_variant"}
REGULATORY_TERMS = {
    "upstream_gene_variant",
    "downstream_gene_variant",
    "regulatory_region_variant",
    "regulatory_region_ablation",
    "TF_binding_site_variant",
    "TFBS_ablation",
    "enhancer_ablation",
    "enhancer_variant",
    "promoter_variant",
    "miRNA",
    "NMD_transcript_variant",
    "non_coding_transcript_exon_variant",
    "intergenic_variant",
}

AD_KEYWORDS = (
    "alzheimer",
    "dementia",
    "cognitive",
    "amyloid",
    "tau ",
    "neurodegener",
)
PD_KEYWORDS = ("parkinson",)
LIPID_KEYWORDS = (
    "cholesterol",
    "ldl",
    "hdl",
    "triglyceride",
    "lipid",
    "apolipoprotein",
    "hyperlipidemia",
    "atherosclerosis",
    "coronary",
    "cetp",
    "statin",
)

ALPHAMISSENSE_PRED_MAP = {
    "B": "benign",
    "A": "ambiguous",
    "P": "pathogenic",
}

app = typer.Typer(add_completion=False, help=__doc__)


@dataclass
class VariantContext:
    rsid: str
    validated_genes: list[str] = field(default_factory=list)


def log_step(msg: str, **kwargs: object) -> None:
    parts = [msg] + [f"{k}={v}" for k, v in kwargs.items()]
    typer.echo(" | ".join(parts))


def load_validated_variants(path: Path) -> tuple[list[VariantContext], pd.DataFrame]:
    if not path.is_file():
        raise FileNotFoundError(f"Input not found: {path}")

    frame = pd.read_csv(path)
    if "canonical_rsid" not in frame.columns:
        raise ValueError(f"Expected column 'canonical_rsid' in {path.name}")

    gene_col = "gene_symbol" if "gene_symbol" in frame.columns else None
    by_rsid: dict[str, list[str]] = {}
    for _, row in frame.iterrows():
        raw = row.get("canonical_rsid")
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            continue
        rsid = str(raw).strip()
        if not RSID_RE.match(rsid):
            continue
        key = rsid.lower()
        genes = by_rsid.setdefault(key, [])
        if gene_col and pd.notna(row.get(gene_col)):
            gene = str(row[gene_col]).strip()
            if gene and gene not in genes:
                genes.append(gene)

    contexts = [
        VariantContext(rsid=rsid, validated_genes=by_rsid[rsid.lower()])
        for rsid in sorted(by_rsid, key=lambda k: by_rsid[k][0] if by_rsid[k] else k)
    ]
    # Preserve original rsID casing from first occurrence in file
    rsid_order: list[str] = []
    seen: set[str] = set()
    for _, row in frame.iterrows():
        raw = row.get("canonical_rsid")
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            continue
        rsid = str(raw).strip()
        if not RSID_RE.match(rsid):
            continue
        key = rsid.lower()
        if key not in seen:
            seen.add(key)
            rsid_order.append(rsid)

    contexts = [
        VariantContext(rsid=rsid, validated_genes=by_rsid[rsid.lower()])
        for rsid in rsid_order
    ]
    return contexts, frame


def cache_path(cache_dir: Path, source: str, rsid: str) -> Path:
    safe = quote(rsid, safe="")
    return cache_dir / source / f"{safe}.json"


def http_get_json(
    session: requests.Session,
    url: str,
    *,
    min_interval_sec: float,
    timeout_sec: float,
    max_retries: int,
    last_request_end: float | None,
) -> tuple[Any | None, float | None, int | None]:
    if last_request_end is not None:
        elapsed = time.monotonic() - last_request_end
        wait = min_interval_sec - elapsed
        if wait > 0:
            time.sleep(wait)

    attempt = 0
    while True:
        attempt += 1
        try:
            response = session.get(url, timeout=timeout_sec)
        except requests.RequestException as exc:
            log_step("request error", url=url, error=str(exc))
            return None, time.monotonic(), None

        if response.status_code == 404:
            return None, time.monotonic(), 404

        if response.status_code in {429, 503}:
            retry_after = response.headers.get("Retry-After")
            sleep_s = min_interval_sec * (2 ** (attempt - 1))
            if retry_after is not None:
                try:
                    sleep_s = max(float(retry_after), sleep_s)
                except ValueError:
                    pass
            if attempt > max_retries:
                log_step("giving up", url=url, status=response.status_code)
                return None, time.monotonic(), response.status_code
            time.sleep(sleep_s)
            continue

        if not response.ok:
            log_step("HTTP error", url=url, status=response.status_code, body=response.text[:200])
            return None, time.monotonic(), response.status_code

        return response.json(), time.monotonic(), response.status_code


def load_or_fetch_json(
    session: requests.Session,
    url: str,
    cache_file: Path,
    *,
    min_interval_sec: float,
    timeout_sec: float,
    max_retries: int,
    last_request_end: float | None,
    legacy_paths: list[Path] | None = None,
) -> tuple[Any | None, bool, float | None]:
    if cache_file.is_file():
        with cache_file.open(encoding="utf-8") as handle:
            return json.load(handle), True, last_request_end

    if legacy_paths:
        for legacy in legacy_paths:
            if legacy.is_file():
                with legacy.open(encoding="utf-8") as handle:
                    payload = json.load(handle)
                cache_file.parent.mkdir(parents=True, exist_ok=True)
                with cache_file.open("w", encoding="utf-8") as handle:
                    json.dump(payload, handle, indent=2)
                return payload, True, last_request_end

    payload, last_request_end, status = http_get_json(
        session,
        url,
        min_interval_sec=min_interval_sec,
        timeout_sec=timeout_sec,
        max_retries=max_retries,
        last_request_end=last_request_end,
    )
    if payload is not None:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        with cache_file.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
    elif status != 404:
        return None, False, last_request_end
    return payload, False, last_request_end


def vep_has_cadd(payload: list[dict[str, Any]]) -> bool:
    if not payload:
        return False
    transcripts = payload[0].get("transcript_consequences") or []
    return any(tc.get("cadd_phred") is not None for tc in transcripts)


def fetch_vep(
    session: requests.Session,
    rsid: str,
    cache_dir: Path,
    last_request_end: float | None,
) -> tuple[list[dict[str, Any]] | None, float | None]:
    url = (
        f"{ENSEMBL_BASE}/vep/human/id/{quote(rsid, safe='')}"
        f"?content-type=application/json&CADD=1"
    )
    cache_file = cache_path(cache_dir, "vep", rsid)
    legacy = LEGACY_VEP_CACHE / f"{rsid}.json"

    if cache_file.is_file():
        with cache_file.open(encoding="utf-8") as handle:
            payload = json.load(handle)
        if isinstance(payload, list) and vep_has_cadd(payload):
            return payload, last_request_end
        cache_file.unlink(missing_ok=True)

    legacy_paths: list[Path] = []
    if legacy.is_file():
        with legacy.open(encoding="utf-8") as handle:
            legacy_payload = json.load(handle)
        if isinstance(legacy_payload, list) and vep_has_cadd(legacy_payload):
            cache_file.parent.mkdir(parents=True, exist_ok=True)
            with cache_file.open("w", encoding="utf-8") as handle:
                json.dump(legacy_payload, handle, indent=2)
            return legacy_payload, last_request_end

    payload, _, last_request_end = load_or_fetch_json(
        session,
        url,
        cache_file,
        min_interval_sec=REQUEST_DELAY_SEC,
        timeout_sec=REQUEST_TIMEOUT_SEC,
        max_retries=MAX_RETRIES,
        last_request_end=last_request_end,
    )
    if payload is None:
        return None, last_request_end
    if isinstance(payload, list):
        return payload, last_request_end
    return [payload], last_request_end


def fetch_myvariant(
    session: requests.Session,
    rsid: str,
    cache_dir: Path,
    last_request_end: float | None,
) -> tuple[dict[str, Any] | list[dict[str, Any]] | None, float | None]:
    fields = "dbnsfp.alphamissense,dbnsfp.alphamissense_pred,cadd.phred,cadd.raw"
    url = f"{MYVARIANT_BASE}/variant/{quote(rsid, safe='')}?fields={fields}"
    cache_file = cache_path(cache_dir, "myvariant", rsid)
    payload, _, last_request_end = load_or_fetch_json(
        session,
        url,
        cache_file,
        min_interval_sec=REQUEST_DELAY_SEC,
        timeout_sec=REQUEST_TIMEOUT_SEC,
        max_retries=MAX_RETRIES,
        last_request_end=last_request_end,
    )
    if payload is None:
        return None, last_request_end
    if isinstance(payload, (dict, list)):
        return payload, last_request_end
    return None, last_request_end


def pick_myvariant_record(
    payload: dict[str, Any] | list[dict[str, Any]] | None,
    ref: str,
    alt: str,
) -> dict[str, Any] | None:
    if payload is None:
        return None
    if isinstance(payload, dict):
        return payload
    if not payload:
        return None

    primary_alt = alt.split(",")[0].strip() if alt else ""
    if primary_alt and ref:
        needle = f"{ref}>{primary_alt}"
        for item in payload:
            if needle in str(item.get("_id", "")):
                return item

    for item in payload:
        am = (item.get("dbnsfp") or {}).get("alphamissense") or {}
        if am.get("score") is not None:
            return item
    return payload[0]


def fetch_gwas_traits_for_association(
    session: requests.Session,
    association_href: str,
    cache_dir: Path,
    last_request_end: float | None,
) -> tuple[list[str], float | None]:
    assoc_id = association_href.rstrip("/").split("/")[-1].split("{")[0]
    cache_file = cache_dir / "gwas_traits" / f"{assoc_id}.json"
    url = f"{GWAS_BASE}/associations/{assoc_id}/efoTraits"
    payload, _, last_request_end = load_or_fetch_json(
        session,
        url,
        cache_file,
        min_interval_sec=REQUEST_DELAY_SEC,
        timeout_sec=REQUEST_TIMEOUT_SEC,
        max_retries=MAX_RETRIES,
        last_request_end=last_request_end,
    )
    if not payload:
        return [], last_request_end
    traits = payload.get("_embedded", {}).get("efoTraits", [])
    names = sorted({str(t.get("trait", "")).strip() for t in traits if t.get("trait")})
    return names, last_request_end


def fetch_gwas_associations(
    session: requests.Session,
    rsid: str,
    cache_dir: Path,
    last_request_end: float | None,
    *,
    max_pages: int = 10,
) -> tuple[list[dict[str, Any]], float | None]:
    cache_file = cache_path(cache_dir, "gwas", rsid)
    if cache_file.is_file():
        with cache_file.open(encoding="utf-8") as handle:
            cached = json.load(handle)
        return cached.get("associations", []), last_request_end

    associations: list[dict[str, Any]] = []
    page = 0
    while page < max_pages:
        url = (
            f"{GWAS_BASE}/singleNucleotidePolymorphisms/{quote(rsid, safe='')}"
            f"/associations?size=100&page={page}"
        )
        payload, last_request_end, status = http_get_json(
            session,
            url,
            min_interval_sec=REQUEST_DELAY_SEC,
            timeout_sec=REQUEST_TIMEOUT_SEC,
            max_retries=MAX_RETRIES,
            last_request_end=last_request_end,
        )
        if payload is None:
            break
        batch = payload.get("_embedded", {}).get("associations", [])
        if not batch:
            break

        for assoc in batch:
            assoc_href = assoc.get("_links", {}).get("self", {}).get("href", "")
            traits, last_request_end = fetch_gwas_traits_for_association(
                session, assoc_href, cache_dir, last_request_end
            )
            p_m = assoc.get("pvalueMantissa")
            p_e = assoc.get("pvalueExponent")
            pvalue = None
            if p_m is not None and p_e is not None:
                pvalue = float(p_m) * (10 ** int(p_e))

            reported_genes: list[str] = []
            for locus in assoc.get("loci", []):
                for gene in locus.get("authorReportedGenes", []):
                    name = gene.get("geneName")
                    if name:
                        reported_genes.append(str(name))

            associations.append(
                {
                    "traits": traits,
                    "pvalue": pvalue,
                    "beta_direction": assoc.get("betaDirection"),
                    "reported_genes": sorted(set(reported_genes)),
                    "association_id": assoc_href.rstrip("/").split("/")[-1].split("{")[0],
                }
            )

        links = payload.get("_links", {})
        if "next" not in links:
            break
        page += 1

    cache_file.parent.mkdir(parents=True, exist_ok=True)
    with cache_file.open("w", encoding="utf-8") as handle:
        json.dump({"rsid": rsid, "associations": associations}, handle, indent=2)
    return associations, last_request_end


def parse_alleles(allele_string: str) -> tuple[str, str]:
    alleles = [a.strip() for a in str(allele_string).split("/") if a.strip()]
    if not alleles:
        return "", ""
    ref = alleles[0]
    alt = alleles[1] if len(alleles) > 1 else ""
    return ref, alt


def pick_canonical_transcript(
    transcripts: list[dict[str, Any]],
    validated_genes: list[str],
) -> dict[str, Any] | None:
    if not transcripts:
        return None

    validated_upper = {g.upper() for g in validated_genes}

    def rank_key(tc: dict[str, Any]) -> tuple[int, int, int, str]:
        gene = str(tc.get("gene_symbol", "")).upper()
        gene_match = 0 if gene in validated_upper else 1
        impact = IMPACT_RANK.get(str(tc.get("impact", "MODIFIER")), 99)
        terms = tc.get("consequence_terms") or []
        has_missense = 0 if "missense_variant" in terms else 1
        return (gene_match, impact, has_missense, str(tc.get("transcript_id", "")))

    return sorted(transcripts, key=rank_key)[0]


def classify_functional(consequence_terms: list[str]) -> str:
    terms = set(consequence_terms)
    if terms & MISSENSE_TERMS:
        return "coding-missense"
    if terms & SYNONYMOUS_TERMS:
        return "coding-synonymous"
    if terms & UTR_TERMS:
        return "UTR"
    if terms & INTRONIC_TERMS:
        return "intronic"
    if terms & REGULATORY_TERMS:
        return "upstream-downstream-regulatory"
    return "other"


def alphagenome_eligible(functional_class: str) -> bool:
    return functional_class in {
        "UTR",
        "intronic",
        "upstream-downstream-regulatory",
        "other",
    }


def extract_alphamissense(
    myvariant_payload: dict[str, Any] | list[dict[str, Any]] | None,
    ref: str,
    alt: str,
) -> tuple[str, str]:
    record = pick_myvariant_record(myvariant_payload, ref, alt)
    if not record:
        return "", ""

    dbnsfp = record.get("dbnsfp") or {}
    am = dbnsfp.get("alphamissense") or {}
    scores = am.get("score")
    preds = am.get("pred") or dbnsfp.get("alphamissense_pred")

    if scores is None:
        return "", ""

    if isinstance(scores, (int, float)):
        score_list = [float(scores)]
    else:
        score_list = [float(s) for s in scores]

    if preds is None:
        pred_list: list[str] = []
    elif isinstance(preds, str):
        pred_list = [preds]
    else:
        pred_list = [str(p) for p in preds]

    idx = 0
    score_val = score_list[idx] if score_list else ""
    score_str = f"{score_val:.4f}" if score_val != "" else ""

    pred_raw = pred_list[idx] if idx < len(pred_list) else (pred_list[0] if pred_list else "")
    pred_class = ALPHAMISSENSE_PRED_MAP.get(str(pred_raw).upper(), str(pred_raw) if pred_raw else "")
    return score_str, pred_class


def extract_cadd(
    canonical_tc: dict[str, Any] | None,
    myvariant_payload: dict[str, Any] | list[dict[str, Any]] | None,
    ref: str,
    alt: str,
) -> str:
    if canonical_tc and canonical_tc.get("cadd_phred") is not None:
        return str(canonical_tc["cadd_phred"])
    record = pick_myvariant_record(myvariant_payload, ref, alt)
    if record:
        cadd = record.get("cadd") or {}
        if cadd.get("phred") is not None:
            return str(cadd["phred"])
    return ""


def trait_matches(trait: str, keywords: tuple[str, ...]) -> bool:
    lower = trait.lower()
    return any(kw in lower for kw in keywords)


def summarize_gwas_traits(associations: list[dict[str, Any]]) -> tuple[str, str, str, str]:
    all_traits: set[str] = set()
    ad_traits: set[str] = set()
    pd_traits: set[str] = set()
    lipid_traits: set[str] = set()

    for assoc in associations:
        for trait in assoc.get("traits", []):
            if not trait:
                continue
            all_traits.add(trait)
            if trait_matches(trait, AD_KEYWORDS):
                ad_traits.add(trait)
            if trait_matches(trait, PD_KEYWORDS):
                pd_traits.add(trait)
            if trait_matches(trait, LIPID_KEYWORDS):
                lipid_traits.add(trait)

    def join(items: set[str]) -> str:
        return "; ".join(sorted(items))

    return join(all_traits), join(ad_traits), join(pd_traits), join(lipid_traits)


def annotate_variant(
    ctx: VariantContext,
    session: requests.Session,
    cache_dir: Path,
    last_request_end: float | None,
    *,
    skip_gwas: bool,
) -> tuple[dict[str, Any], float | None]:
    rsid = ctx.rsid
    vep_payload, last_request_end = fetch_vep(session, rsid, cache_dir, last_request_end)
    myvariant_payload, last_request_end = fetch_myvariant(session, rsid, cache_dir, last_request_end)

    gwas_associations: list[dict[str, Any]] = []
    if not skip_gwas:
        gwas_associations, last_request_end = fetch_gwas_associations(
            session, rsid, cache_dir, last_request_end
        )

    if not vep_payload:
        all_traits, ad_traits, pd_traits, lipid_traits = summarize_gwas_traits(gwas_associations)
        return {
            "rsID": rsid,
            "genome_build": GENOME_BUILD,
            "validated_gene_symbols": ";".join(ctx.validated_genes),
            "chromosome_GRCh38": "NA",
            "position_GRCh38": "NA",
            "ref_allele": "NA",
            "alt_allele": "NA",
            "most_severe_consequence": "NA",
            "canonical_gene": "NA",
            "canonical_transcript": "NA",
            "canonical_consequence_terms": "NA",
            "functional_class": "NA",
            "alphagenome_eligible": False,
            "SIFT": "NA",
            "PolyPhen": "NA",
            "CADD_phred": "NA",
            "AlphaMissense_score": "NA",
            "AlphaMissense_class": "NA",
            "no_vep_record": True,
            "gwas_all_traits": all_traits or "NA",
            "gwas_ad_traits": ad_traits or "NA",
            "gwas_pd_traits": pd_traits or "NA",
            "gwas_lipid_traits": lipid_traits or "NA",
            "gwas_association_count": len(gwas_associations),
        }, last_request_end

    variant = vep_payload[0]
    transcripts = variant.get("transcript_consequences") or []
    if not isinstance(transcripts, list):
        transcripts = []

    ref, alt = parse_alleles(str(variant.get("allele_string", "")))
    canonical = pick_canonical_transcript(transcripts, ctx.validated_genes)
    canonical_terms = list(canonical.get("consequence_terms") or []) if canonical else []
    functional_class = classify_functional(canonical_terms) if canonical_terms else "other"

    sift = str(canonical.get("sift_prediction") or canonical.get("sift") or "") if canonical else ""
    polyphen = str(canonical.get("polyphen_prediction") or canonical.get("polyphen") or "") if canonical else ""
    cadd = extract_cadd(canonical, myvariant_payload, ref, alt)
    am_score, am_class = extract_alphamissense(myvariant_payload, ref, alt)

    if functional_class != "coding-missense":
        am_score, am_class = "", ""

    all_traits, ad_traits, pd_traits, lipid_traits = summarize_gwas_traits(gwas_associations)

    return {
        "rsID": rsid,
        "genome_build": GENOME_BUILD,
        "validated_gene_symbols": ";".join(ctx.validated_genes),
        "chromosome_GRCh38": str(variant.get("seq_region_name", "")),
        "position_GRCh38": str(variant.get("start", "")),
        "ref_allele": ref,
        "alt_allele": alt,
        "most_severe_consequence": str(variant.get("most_severe_consequence", "")),
        "canonical_gene": str(canonical.get("gene_symbol", "")) if canonical else "",
        "canonical_transcript": str(canonical.get("transcript_id", "")) if canonical else "",
        "canonical_consequence_terms": ",".join(canonical_terms),
        "functional_class": functional_class,
        "alphagenome_eligible": alphagenome_eligible(functional_class),
        "SIFT": sift or "NA",
        "PolyPhen": polyphen or "NA",
        "CADD_phred": cadd or "NA",
        "AlphaMissense_score": am_score or "NA",
        "AlphaMissense_class": am_class or "NA",
        "no_vep_record": False,
        "gwas_all_traits": all_traits or "NA",
        "gwas_ad_traits": ad_traits or "NA",
        "gwas_pd_traits": pd_traits or "NA",
        "gwas_lipid_traits": lipid_traits or "NA",
        "gwas_association_count": len(gwas_associations),
    }, last_request_end


def write_summary_md(
    output_path: Path,
    annotation_df: pd.DataFrame,
    input_df: pd.DataFrame,
    unresolved_count: int,
) -> None:
    total_unique = len(annotation_df)
    class_counts = annotation_df["functional_class"].value_counts().to_dict()
    coding_missense = int(class_counts.get("coding-missense", 0))
    coding_synonymous = int(class_counts.get("coding-synonymous", 0))
    coding_total = coding_missense + coding_synonymous
    noncoding_total = total_unique - coding_total - int(class_counts.get("NA", 0))
    no_record = int(annotation_df["no_vep_record"].sum())
    alphagenome_ready = annotation_df[annotation_df["alphagenome_eligible"]].copy()

    lines = [
        "# Coding vs non-coding variant summary",
        "",
        f"**Genome build:** {GENOME_BUILD}",
        f"**Input:** `results/snps_validated.csv` ({len(input_df)} rows; "
        f"{total_unique} unique resolved rsIDs; {unresolved_count} rows without rsID)",
        "",
        "## Functional class counts",
        "",
        "| Class | Count |",
        "|---|---:|",
    ]
    for cls in [
        "coding-missense",
        "coding-synonymous",
        "UTR",
        "intronic",
        "upstream-downstream-regulatory",
        "other",
        "NA",
    ]:
        count = int(class_counts.get(cls, 0))
        if count or cls != "NA":
            lines.append(f"| {cls} | {count} |")

    lines.extend(
        [
            "",
            "## Manuscript interpretability check",
            "",
            f"- **Coding variants (missense + synonymous):** {coding_total} "
            f"({coding_missense} missense, {coding_synonymous} synonymous)",
            f"- **Non-coding / regulatory variants:** {noncoding_total}",
            f"- **No Ensembl VEP record:** {no_record}",
            "",
            "AlphaGenome predicts **regulatory / expression-level** effects from sequence context. "
            "Its \"no expression change\" result is therefore **only interpretable for the "
            f"{len(alphagenome_ready)} non-coding variants** listed below — **not** for the "
            f"{coding_total} coding variants.",
            "",
            "**Coding missense variants** ("
            f"{coding_missense}) should be interpreted with "
            "**AlphaMissense** (benign / ambiguous / pathogenic) plus SIFT, PolyPhen, and CADD — "
            "these capture qualitative protein effects that AlphaGenome does not model.",
            "",
            f"**Coding synonymous variants** ({coding_synonymous}) change codon usage without amino-acid "
            "substitution; they are excluded from the AlphaGenome-ready set for the same reason as "
            "missense variants (AlphaGenome is not the primary evidence for protein-level claims).",
            "",
            "## AlphaGenome-ready variants (non-coding subset)",
            "",
            "Format matches the existing AlphaGenome pipeline (`rsID`, GRCh38 chromosome, position, ref, alt).",
            "",
            "| rsID | chr | GRCh38 pos | ref | alt | functional_class | validated_gene(s) |",
            "|---|---|---:|---|---|---|---|",
        ]
    )

    for _, row in alphagenome_ready.sort_values("rsID").iterrows():
        lines.append(
            f"| {row['rsID']} | {row['chromosome_GRCh38']} | {row['position_GRCh38']} | "
            f"{row['ref_allele']} | {row['alt_allele']} | {row['functional_class']} | "
            f"{row['validated_gene_symbols']} |"
        )

    lines.extend(
        [
            "",
            "## Coding missense — AlphaMissense summary",
            "",
            "| rsID | gene | AlphaMissense class | score | SIFT | PolyPhen | CADD |",
            "|---|---|---|---:|---|---|---:|",
        ]
    )
    missense_df = annotation_df[annotation_df["functional_class"] == "coding-missense"].sort_values("rsID")
    for _, row in missense_df.iterrows():
        lines.append(
            f"| {row['rsID']} | {row['canonical_gene']} | {row['AlphaMissense_class']} | "
            f"{row['AlphaMissense_score']} | {row['SIFT']} | {row['PolyPhen']} | {row['CADD_phred']} |"
        )

    lines.extend(
        [
            "",
            "## Unresolved input rows (no rsID)",
            "",
        ]
    )
    unresolved = input_df[
        input_df["canonical_rsid"].isna()
        | ~input_df["canonical_rsid"].astype(str).str.match(RSID_RE, na=False)
    ]
    if unresolved.empty:
        lines.append("_None._")
    else:
        lines.append("| gene_symbol | snp_identifier_raw | resolution_notes |")
        lines.append("|---|---|---|")
        for _, row in unresolved.iterrows():
            lines.append(
                f"| {row.get('gene_symbol', '')} | {row.get('snp_identifier_raw', '')} | "
                f"{row.get('resolution_notes', '')} |"
            )

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@app.callback(invoke_without_command=True)
def main(
    input_csv: Path = typer.Option(DEFAULT_INPUT, "--input", help="Validated SNP table."),
    output_dir: Path = typer.Option(DEFAULT_OUTPUT_DIR, "--output-dir"),
    cache_dir: Path = typer.Option(DEFAULT_CACHE_DIR, "--cache-dir"),
    skip_gwas: bool = typer.Option(False, "--skip-gwas", help="Skip GWAS Catalog lookups."),
) -> None:
    """Annotate unique validated rsIDs and write CSV + summary markdown."""
    contexts, input_df = load_validated_variants(input_csv)
    if not contexts:
        raise typer.Exit("No resolved rsIDs found in input.")

    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update(
        {
            "Content-Type": "application/json",
            "User-Agent": "rogen-aging-variant-annotation/1.0 (academic research)",
        }
    )

    rows: list[dict[str, Any]] = []
    last_request_end: float | None = None

    for index, ctx in enumerate(contexts, start=1):
        log_step(f"[{index}/{len(contexts)}]", rsid=ctx.rsid)
        row, last_request_end = annotate_variant(
            ctx,
            session,
            cache_dir,
            last_request_end,
            skip_gwas=skip_gwas,
        )
        rows.append(row)

    session.close()

    annotation_df = pd.DataFrame(rows)
    csv_path = output_dir / "variant_functional_annotation.csv"
    summary_path = output_dir / "coding_vs_noncoding_summary.md"
    annotation_df.to_csv(csv_path, index=False)

    unresolved_count = int(
        input_df["canonical_rsid"].isna().sum()
        + (
            input_df["canonical_rsid"].notna()
            & ~input_df["canonical_rsid"].astype(str).str.match(RSID_RE, na=False)
        ).sum()
    )
    write_summary_md(summary_path, annotation_df, input_df, unresolved_count)

    log_step("done", annotated=len(annotation_df), csv=str(csv_path), summary=str(summary_path))


if __name__ == "__main__":
    app()
