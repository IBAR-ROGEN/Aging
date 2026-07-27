#!/usr/bin/env python3
"""Cluster ∩ LongevityMap overlap enrichment analysis (GRCh38/hg38 gene symbols).

Tests whether AD/PD cluster genes overlap LongevityMap longevity-associated genes
more than expected under documented universe definitions.
"""

from __future__ import annotations

import json
import re
import sqlite3
import time
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from urllib.parse import quote

import pandas as pd
import requests
import typer
from scipy.stats import fisher_exact, hypergeom

GENOME_BUILD = "GRCh38/hg38"
REPO_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_CLUSTER_TABLE = REPO_ROOT / "data" / "Supplementary Table 3.xlsx"
DEFAULT_LONGEVITY_DB = REPO_ROOT / "data" / "longevitymap.sqlite"
DEFAULT_SNPS_VALIDATED = REPO_ROOT / "results" / "snps_validated.csv"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results"
FALLBACK_SNPS_VALIDATED = REPO_ROOT / "analysis" / "results" / "snps_validated.csv"
FALLBACK_OUTPUT_DIR = REPO_ROOT / "analysis" / "results"
DEFAULT_GPL_FILE = REPO_ROOT / "data" / "ad_pd_gpl_ids.txt"
DEFAULT_CACHE_DIR = DEFAULT_OUTPUT_DIR / "cache"

LONGEVITY_ZIP_URL = "https://www.genomics.senescence.info/longevity/longevity_genes.zip"
GENCODE_GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/"
    "gencode.v46.primary_assembly.annotation.gtf.gz"
)

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

EXPECTED_OVERLAP = 41
REQUEST_DELAY_SEC = 0.34

app = typer.Typer(add_completion=False, help=__doc__)


@app.command("build-sqlite")
def build_sqlite_cmd(
    longevity_db: Path = typer.Option(DEFAULT_LONGEVITY_DB, "--longevity-db"),
    cache_dir: Path = typer.Option(REPO_ROOT / "results" / "cache", "--cache-dir"),
) -> None:
    """Build data/longevitymap.sqlite from the official HAGR LongevityMap CSV."""
    build_longevity_sqlite(longevity_db, cache_dir)


def log_step(msg: str, **kwargs: object) -> None:
    parts = [msg] + [f"{k}={v}" for k, v in kwargs.items()]
    typer.echo(" | ".join(parts))


def normalize_symbol(symbol: str) -> str:
    return str(symbol).strip().upper()


def split_gene_field(raw: str) -> list[str]:
    genes: list[str] = []
    for part in re.split(r"[,;/]", str(raw)):
        token = part.strip()
        if token:
            genes.append(normalize_symbol(token))
    return genes


def read_cache_json(cache_dir: Path, key: str) -> dict | list | None:
    path = cache_dir / f"{key}.json"
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def write_cache_json(cache_dir: Path, key: str, payload: dict | list) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{key}.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def download_longevity_csv(cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    zip_path = cache_dir / "longevity_genes.zip"
    csv_path = cache_dir / "longevity.csv"
    if csv_path.is_file():
        return csv_path
    typer.echo(f"Downloading LongevityMap CSV from {LONGEVITY_ZIP_URL}")
    resp = requests.get(LONGEVITY_ZIP_URL, timeout=120)
    resp.raise_for_status()
    zip_path.write_bytes(resp.content)
    with zipfile.ZipFile(zip_path) as zf:
        with zf.open("longevity.csv") as src, csv_path.open("wb") as dst:
            dst.write(src.read())
    log_step("Downloaded LongevityMap CSV", rows_out=sum(1 for _ in csv_path.open()) - 1)
    return csv_path


def build_longevity_sqlite(db_path: Path, cache_dir: Path) -> None:
    csv_path = download_longevity_csv(cache_dir)
    raw = pd.read_csv(csv_path)
    raw.columns = [c.strip() for c in raw.columns]

    variant_rows: list[dict[str, object]] = []
    gene_symbols: set[str] = set()

    for _, row in raw.iterrows():
        association = str(row.get("Association", "")).strip().lower()
        if association == "non-significant":
            association = "non-significant"
        genes = split_gene_field(row.get("Gene(s)", ""))
        for gene in genes:
            gene_symbols.add(gene)
        variant_rows.append(
            {
                "id": row.get("id"),
                "Association": association,
                "Population": row.get("Population"),
                "Variant(s)": row.get("Variant(s)"),
                "Gene(s)": row.get("Gene(s)"),
                "PubMed": row.get("PubMed"),
            }
        )

    variant_df = pd.DataFrame(variant_rows)
    gene_df = pd.DataFrame({"Gene Symbol": sorted(gene_symbols)})

    db_path.parent.mkdir(parents=True, exist_ok=True)
    if db_path.is_file():
        db_path.unlink()
    conn = sqlite3.connect(db_path)
    variant_df.to_sql("variant", conn, index=False, if_exists="replace")
    gene_df.to_sql("gene", conn, index=False, if_exists="replace")
    conn.close()
    log_step(
        "Built longevitymap.sqlite",
        variants=len(variant_df),
        unique_genes=len(gene_df),
        path=str(db_path),
    )



def longevity_gene_sets(db_path: Path) -> tuple[set[str], set[str]]:
    conn = sqlite3.connect(db_path)
    tables = {
        row[0]
        for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
    }

    if "gene" in tables and "variant" in tables:
        variant_df = pd.read_sql_query(
            """
            SELECT LOWER(v.association) AS association, g.symbol AS symbol
            FROM variant v
            JOIN gene g ON v.gene_id = g.id
            WHERE g.symbol IS NOT NULL AND TRIM(g.symbol) != ''
            """,
            conn,
        )
    else:
        variant_df = pd.read_sql_query("SELECT Association, [Gene(s)] AS genes FROM variant", conn)
        variant_df = variant_df.rename(columns={"Association": "association", "genes": "symbol"})
        variant_df["association"] = variant_df["association"].astype(str).str.lower()
        variant_df = variant_df.assign(
            symbol=variant_df["symbol"].map(lambda x: split_gene_field(x)[0] if split_gene_field(x) else "")
        )
    conn.close()

    sig_genes: set[str] = set()
    all_genes: set[str] = set()
    for _, row in variant_df.iterrows():
        symbol = normalize_symbol(row["symbol"])
        if not symbol or symbol == "NA":
            continue
        all_genes.add(symbol)
        if str(row["association"]).strip().lower() == "significant":
            sig_genes.add(symbol)
    log_step(
        "Loaded LongevityMap gene sets from sqlite",
        significant_genes=len(sig_genes),
        all_genes=len(all_genes),
    )
    return sig_genes, all_genes


def fetch_entrez_symbols(entrez_ids: Iterable[int], cache_dir: Path) -> dict[int, str]:
    mapping: dict[int, str] = {}
    pending: list[int] = []
    for eid in sorted(set(int(x) for x in entrez_ids)):
        cached = read_cache_json(cache_dir, f"entrez_symbol_{eid}")
        if cached and "symbol" in cached:
            mapping[eid] = cached["symbol"]
        else:
            pending.append(eid)

    batch_size = 200
    for start in range(0, len(pending), batch_size):
        batch = pending[start : start + batch_size]
        id_str = ",".join(str(eid) for eid in batch)
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=gene&id={id_str}&retmode=json"
        )
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        result = resp.json().get("result", {})
        uids = result.get("uids", [])
        for eid in batch:
            entry = result.get(str(eid), {})
            symbol = normalize_symbol(entry.get("name", "NA"))
            mapping[eid] = symbol
            write_cache_json(cache_dir, f"entrez_symbol_{eid}", {"entrez": eid, "symbol": symbol})
        if len(uids) != len(batch):
            typer.echo(f"WARNING: esummary batch returned {len(uids)} of {len(batch)} IDs")
        time.sleep(REQUEST_DELAY_SEC)

    log_step("Mapped ENTREZ IDs to symbols", mapped=len(mapping), pending=len(pending))
    return mapping


def load_cluster_sets(cluster_path: Path, cache_dir: Path) -> dict[str, set[str]]:
    if not cluster_path.is_file():
        raise FileNotFoundError(
            f"Cluster table not found: {cluster_path}. "
            "Place Supplementary Table 3.xlsx in data/."
        )

    xls = pd.ExcelFile(cluster_path, engine="openpyxl")
    entrez_by_sheet: dict[str, list[int]] = {}

    for key, sheet in CLUSTER_SHEETS.items():
        if sheet not in xls.sheet_names:
            raise ValueError(f"Missing sheet {sheet!r} in {cluster_path}")
        df = pd.read_excel(xls, sheet_name=sheet)
        entrez_col = next((c for c in df.columns if str(c).upper() == "ENTREZ"), None)
        if entrez_col is None:
            symbol_col = next((c for c in df.columns if "symbol" in str(c).lower()), None)
            if symbol_col is None:
                raise ValueError(f"Sheet {sheet!r} lacks ENTREZ or gene symbol column")
        else:
            entrez_by_sheet[key] = [int(x) for x in df[entrez_col].dropna()]

    all_entrez = {eid for ids in entrez_by_sheet.values() for eid in ids}
    symbol_map = fetch_entrez_symbols(all_entrez, cache_dir) if all_entrez else {}

    cluster_sets: dict[str, set[str]] = {}
    for key, sheet in CLUSTER_SHEETS.items():
        df = pd.read_excel(xls, sheet_name=sheet)
        symbol_col = next((c for c in df.columns if "symbol" in str(c).lower()), None)
        symbols: set[str] = set()
        if symbol_col is not None:
            symbols = {normalize_symbol(x) for x in df[symbol_col].dropna() if str(x).strip()}
        if key in entrez_by_sheet:
            mapped = {
                symbol_map[eid]
                for eid in entrez_by_sheet[key]
                if eid in symbol_map and symbol_map[eid] != "NA"
            }
            symbols |= mapped
        cluster_sets[key] = symbols
        stated = CLUSTER_STATED_SIZES[key]
        log_step(
            f"Loaded cluster {key}",
            rows_in=len(df),
            unique_genes=len(symbols),
            stated=stated,
        )
        if len(symbols) != stated:
            typer.echo(
                f"WARNING: cluster {key} size {len(symbols)} != stated {stated} ({GENOME_BUILD})"
            )

    derived = {
        "combined": set().union(*cluster_sets.values()),
        "AD": cluster_sets["AD-up"] | cluster_sets["AD-down"],
        "PD": cluster_sets["PD-up"] | cluster_sets["PD-down"],
        **cluster_sets,
    }
    log_step(
        "Cluster set summary",
        combined=len(derived["combined"]),
        AD=len(derived["AD"]),
        PD=len(derived["PD"]),
    )
    return derived


def fetch_protein_coding_genes(cache_dir: Path) -> set[str]:
    cached = read_cache_json(cache_dir, "gencode_v46_protein_coding_symbols")
    if cached:
        symbols = {normalize_symbol(x) for x in cached}
        log_step("Loaded protein-coding universe (cache)", genes=len(symbols))
        return symbols

    typer.echo(f"Downloading GENCODE GTF from {GENCODE_GTF_URL}")
    resp = requests.get(GENCODE_GTF_URL, timeout=300)
    resp.raise_for_status()
    import gzip
    from io import BytesIO

    symbols: set[str] = set()
    with gzip.GzipFile(fileobj=BytesIO(resp.content)) as gz:
        for line in gz:
            text = line.decode("utf-8").strip()
            if not text or text.startswith("#") or "\tgene\t" not in text:
                continue
            if "gene_type \"protein_coding\"" not in text:
                continue
            match = re.search(r'gene_name "([^"]+)"', text)
            if match:
                symbols.add(normalize_symbol(match.group(1)))
    write_cache_json(cache_dir, "gencode_v46_protein_coding_symbols", sorted(symbols))
    log_step("Built protein-coding universe", genes=len(symbols))
    return symbols


def fetch_platform_genes(gpl_file: Path, cache_dir: Path) -> set[str]:
    if not gpl_file.is_file():
        return set()
    gpl_ids = [
        line.strip()
        for line in gpl_file.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]
    if not gpl_ids:
        return set()

    import GEOparse

    symbols: set[str] = set()
    for gpl in gpl_ids:
        cache_key = f"gpl_genes_{gpl}"
        cached = read_cache_json(cache_dir, cache_key)
        if cached:
            gpl_symbols = {normalize_symbol(x) for x in cached}
        else:
            gse_gpl = GEOparse.get_GEO(geo=gpl, destdir=str(cache_dir / "geo_gpl"))
            gpl_symbols = set()
            table = getattr(gse_gpl, "table", None)
            if table is not None and "Gene Symbol" in table.columns:
                for val in table["Gene Symbol"].dropna():
                    for token in re.split(r"[/]{3}|,", str(val)):
                        token = token.strip()
                        if token and token.lower() != "na":
                            gpl_symbols.add(normalize_symbol(token))
            write_cache_json(cache_dir, cache_key, sorted(gpl_symbols))
            time.sleep(REQUEST_DELAY_SEC)
        symbols |= gpl_symbols
        log_step(f"Loaded platform {gpl}", genes=len(gpl_symbols))
    log_step("Microarray platform union universe", genes=len(symbols), platforms=len(gpl_ids))
    return symbols


@dataclass
class EnrichmentResult:
    cluster_label: str
    b_label: str
    universe_label: str
    universe_size: int
    cluster_in_universe: int
    longevity_in_universe: int
    overlap: int
    expected: float
    fold_enrichment: float
    odds_ratio: float
    fisher_p: float
    hypergeom_p: float
    bh_q_fisher: float | None = None
    bh_q_hypergeom: float | None = None

    def contingency_markdown(self) -> str:
        a = self.cluster_in_universe
        b = self.overlap
        c = a - b
        d = self.longevity_in_universe - b
        e = self.universe_size - self.cluster_in_universe - self.longevity_in_universe + b
        return (
            f"| | In LongevityMap set | Not in LongevityMap set | Total |\n"
            f"|---|---:|---:|---:|\n"
            f"| In cluster set | {b} | {c} | {a} |\n"
            f"| Not in cluster set | {d} | {e} | {self.universe_size - a} |\n"
            f"| Total | {self.longevity_in_universe} | {self.universe_size - self.longevity_in_universe} | {self.universe_size} |"
        )


def enrich_test(
    cluster: set[str],
    longevity: set[str],
    universe: set[str],
    cluster_label: str,
    b_label: str,
    universe_label: str,
) -> EnrichmentResult:
    uni = set(universe)
    a_set = cluster & uni
    b_set = longevity & uni
    overlap = a_set & b_set

    n = len(uni)
    a = len(a_set)
    b = len(b_set)
    k = len(overlap)

    expected = (a * b / n) if n else float("nan")
    fold = (k / a) / (b / n) if a and b and n else float("nan")

    not_a = n - a
    b_not_k = b - k
    not_b_in_universe = n - b
    a_not_k = a - k
    table = [[k, a_not_k], [b_not_k, not_b_in_universe - a_not_k]]
    _, fisher_p = fisher_exact(table, alternative="greater")

    # M=pop size, n=success states, N=draws, k=observed
    hyper_p = float(hypergeom.sf(k - 1, n, b, a)) if n and a and b else float("nan")

    odds = (k / a_not_k) / (b_not_k / (not_b_in_universe - a_not_k)) if a_not_k and b_not_k else float("inf")

    return EnrichmentResult(
        cluster_label=cluster_label,
        b_label=b_label,
        universe_label=universe_label,
        universe_size=n,
        cluster_in_universe=a,
        longevity_in_universe=b,
        overlap=k,
        expected=expected,
        fold_enrichment=fold,
        odds_ratio=odds,
        fisher_p=fisher_p,
        hypergeom_p=hyper_p,
    )


def bh_correction(p_values: list[float]) -> list[float]:
    m = len(p_values)
    if m == 0:
        return []
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    q_vals = [1.0] * m
    prev = 1.0
    for rank, (idx, p) in enumerate(reversed(indexed), start=1):
        i = m - rank
        q = min(prev, p * m / (i + 1))
        q_vals[idx] = q
        prev = q
    return q_vals


def resolve_snps_validated(path: Path) -> Path:
    if path.is_file():
        return path
    if FALLBACK_SNPS_VALIDATED.is_file():
        return FALLBACK_SNPS_VALIDATED
    return path


def resolve_output_dir(path: Path) -> Path:
    if path == DEFAULT_OUTPUT_DIR and not path.is_dir() and FALLBACK_OUTPUT_DIR.is_dir():
        return FALLBACK_OUTPUT_DIR
    return path


def load_snps_validated_genes(path: Path) -> set[str]:
    resolved = resolve_snps_validated(path)
    if not resolved.is_file():
        return set()
    df = pd.read_csv(resolved)
    assoc_col = next(
        (c for c in df.columns if str(c).lower().replace(" ", "_") in {"snp_association", "snpassociation"}),
        None,
    )
    if assoc_col is None:
        assoc_col = "SNP Association" if "SNP Association" in df.columns else df.columns[0]
    sig = df[df[assoc_col].astype(str).str.strip().str.lower() == "significant"]
    col = "gene_symbol" if "gene_symbol" in sig.columns else "Gene Symbol"
    genes = {normalize_symbol(x) for x in sig[col].dropna()}
    log_step("Loaded snps_validated.csv significant genes", genes=len(genes), path=str(resolved))
    return genes


def render_report(
    results: list[EnrichmentResult],
    overlap_check: dict[str, object],
    output_path: Path,
    primary_universe: str,
) -> None:
    sig_primary = [
        r
        for r in results
        if r.universe_label == primary_universe and r.b_label == "significant"
    ]
    any_sig = any(r.bh_q_fisher is not None and r.bh_q_fisher < 0.05 for r in sig_primary)

    lines = [
        f"# Cluster ∩ LongevityMap overlap enrichment ({GENOME_BUILD})",
        "",
        "Gene symbols are matched on HGNC uppercase symbols. LongevityMap gene sets are",
        "derived from `data/longevitymap.sqlite` (HAGR LongevityMap build 3 CSV).",
        "",
        "## Overlap sanity check (|A ∩ B| vs 41)",
        "",
        f"- Expected overlap (manuscript): **{EXPECTED_OVERLAP}**",
        f"- Computed |A ∩ B| (combined clusters ∩ LongevityMap significant): **{overlap_check.get('computed_overlap', 'NA')}**",
        f"- snps_validated.csv significant unique genes: **{overlap_check.get('snps_validated_genes', 'NA')}**",
        f"- Verdict: **{overlap_check.get('verdict', 'NA')}**",
        "",
    ]
    if overlap_check.get("missing_in_overlap"):
        lines.append(
            f"- Genes in snps_validated but not in computed overlap: {overlap_check['missing_in_overlap']}"
        )
    if overlap_check.get("extra_in_overlap"):
        lines.append(
            f"- Genes in computed overlap but not in snps_validated: {overlap_check['extra_in_overlap']}"
        )
    lines.extend(
        [
            "",
            f"Primary universe for interpretation: **{primary_universe}**",
            "",
            "## Summary",
            "",
        ]
    )

    if not any_sig:
        lines.append(
            "Under the primary honest background universe, **no cluster set shows significant "
            "overlap with LongevityMap genes after Benjamini–Hochberg correction** "
            "(Fisher exact, one-sided enrichment)."
        )
    else:
        lines.append(
            "Some cluster sets show BH-significant enrichment under the primary universe; "
            "see tables below."
        )

    lines.extend(["", "## All tests", ""])
    for result in results:
        sig = (
            "significant"
            if result.bh_q_fisher is not None and result.bh_q_fisher < 0.05
            else "not significant"
        )
        lines.extend(
            [
                f"### {result.cluster_label} × {result.b_label} — universe {result.universe_label}",
                "",
                result.contingency_markdown(),
                "",
                f"- Overlap: **{result.overlap}**",
                f"- Expected by chance: **{result.expected:.3f}**",
                f"- Fold enrichment: **{result.fold_enrichment:.3f}**",
                f"- Odds ratio: **{result.odds_ratio:.3f}**",
                f"- Fisher p: **{result.fisher_p:.4g}** (BH q={result.bh_q_fisher:.4g})",
                f"- Hypergeometric p: **{result.hypergeom_p:.4g}** (BH q={result.bh_q_hypergeom:.4g})",
                f"- Verdict (BH q<0.05, Fisher): **{sig}**",
                "",
            ]
        )

    lines.extend(["", "## Plain-English interpretation", ""])
    lines.append(
        f"We asked whether genes dysregulated in AD/PD brain meta-analysis clusters overlap "
        f"longevity-associated genes from LongevityMap more than expected. "
        f"The combined cluster list shares **{overlap_check.get('computed_overlap', 'NA')}** genes "
        f"with LongevityMap significant entries (manuscript cites **{EXPECTED_OVERLAP}**). "
    )
    if primary_universe == "microarray_platform_union":
        lines.append(
            "The primary background universe is the union of microarray platform gene annotations, "
            "which reflects genes that could have been detected on the arrays used in the meta-analysis."
        )
    else:
        lines.append(
            "The primary background universe is the set of all genes appearing in any cluster DE list, "
            "which reflects genes tested in the meta-analysis."
        )
    if not any_sig:
        lines.append(
            " After correcting for multiple testing, the overlap is **not statistically significant** "
            "under that honest background — the observed overlap is consistent with chance sampling."
        )
    else:
        lines.append(
            " After BH correction, at least one cluster subset remains enriched; inspect per-table q-values."
        )
    lines.append("")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    log_step("Wrote overlap enrichment report", path=str(output_path))


@app.callback(invoke_without_command=True)
def run_enrichment(
    ctx: typer.Context,
    cluster_table: Path = typer.Option(DEFAULT_CLUSTER_TABLE, "--cluster-table"),
    longevity_db: Path = typer.Option(DEFAULT_LONGEVITY_DB, "--longevity-db"),
    snps_validated: Path = typer.Option(DEFAULT_SNPS_VALIDATED, "--snps-validated"),
    gpl_file: Path = typer.Option(DEFAULT_GPL_FILE, "--gpl-file"),
    output_dir: Path = typer.Option(DEFAULT_OUTPUT_DIR, "--output-dir"),
    cache_dir: Path = typer.Option(DEFAULT_OUTPUT_DIR / "cache", "--cache-dir"),
) -> None:
    """Run cluster ∩ LongevityMap overlap enrichment."""
    if ctx.invoked_subcommand is not None:
        return

    typer.echo(f"Genome build: {GENOME_BUILD}")

    snps_validated = resolve_snps_validated(snps_validated)
    output_dir = resolve_output_dir(output_dir)
    cache_dir = output_dir / "cache"

    if not longevity_db.is_file():
        build_longevity_sqlite(longevity_db, cache_dir)

    b_sig, b_all = longevity_gene_sets(longevity_db)
    cluster_sets = load_cluster_sets(cluster_table, cache_dir)

    universe_a = fetch_protein_coding_genes(cache_dir)
    universe_b = fetch_platform_genes(gpl_file, cache_dir)
    universe_c = cluster_sets["combined"]
    universes = {
        "all_protein_coding": universe_a,
        "microarray_platform_union": universe_b,
        "meta_analysis_de_tested": universe_c,
    }
    primary_universe = "microarray_platform_union" if universe_b else "meta_analysis_de_tested"

    cluster_labels = ["combined", "AD", "PD", "AD-up", "AD-down", "PD-up", "PD-down"]
    b_sets = {"significant": b_sig, "all_longevitymap": b_all}

    results: list[EnrichmentResult] = []
    for cl in cluster_labels:
        for b_name, b_genes in b_sets.items():
            for u_name, u_genes in universes.items():
                if u_name == "microarray_platform_union" and not u_genes:
                    continue
                res = enrich_test(cluster_sets[cl], b_genes, u_genes, cl, b_name, u_name)
                results.append(res)

    fisher_ps = [r.fisher_p for r in results]
    hyper_ps = [r.hypergeom_p for r in results]
    fisher_qs = bh_correction(fisher_ps)
    hyper_qs = bh_correction(hyper_ps)
    for r, fq, hq in zip(results, fisher_qs, hyper_qs, strict=True):
        r.bh_q_fisher = fq
        r.bh_q_hypergeom = hq

    computed_overlap = len(cluster_sets["combined"] & b_sig)
    validated_genes = load_snps_validated_genes(snps_validated)
    overlap_genes = cluster_sets["combined"] & b_sig
    verdict = "MATCH" if computed_overlap == EXPECTED_OVERLAP else "MISMATCH"
    overlap_check = {
        "computed_overlap": computed_overlap,
        "snps_validated_genes": len(validated_genes),
        "verdict": verdict,
        "missing_in_overlap": sorted(validated_genes - overlap_genes) if validated_genes else [],
        "extra_in_overlap": sorted(overlap_genes - validated_genes) if validated_genes else [],
    }
    log_step(
        "Overlap check",
        computed=computed_overlap,
        expected=EXPECTED_OVERLAP,
        verdict=verdict,
    )

    stats_csv = output_dir / "overlap_enrichment_stats.csv"
    pd.DataFrame([r.__dict__ for r in results]).to_csv(stats_csv, index=False)
    render_report(results, overlap_check, output_dir / "overlap_enrichment.md", primary_universe)
    typer.echo(f"Done. Report: {output_dir / 'overlap_enrichment.md'}")


if __name__ == "__main__":
    app()
