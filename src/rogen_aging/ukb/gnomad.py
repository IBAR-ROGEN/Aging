"""Activity 2.1.8.1 — Compare 1KG LA-SNP allele frequencies to gnomAD v4 (NFE).

Reads ``analysis/la_snp_1kg_frequencies.csv`` (from ``rogen-ukb-manifest extract``),
fetches matching gnomAD v4 non-Finnish European (``nfe``) allele frequencies via the public GraphQL API, and writes a comparison table plus a
1KG-vs-gnomAD scatter plot.

gnomAD variant lookup uses GRCh38 ``chrom-pos-ref-alt`` IDs built from the 1KG
table (or a single-base region query matched by rsID when alleles are missing).
Responses are cached under ``data/geo/`` so reruns are offline-friendly.

Example:
    uv run rogen-compare-af-gnomad \\
        --input analysis/la_snp_1kg_frequencies.csv \\
        --output analysis/la_snp_af_1kg_vs_gnomad.csv \\
        --scatter figures/af_1kg_vs_gnomad_scatter.png

    uv run rogen-compare-af-gnomad summarize \\
        --input analysis/la_snp_af_1kg_vs_gnomad.csv \\
        --output analysis/af_comparison_summary.md
"""

from __future__ import annotations

import json
import logging
import sys
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Literal, cast

import click
import matplotlib.pyplot as plt
import pandas as pd
import requests
import typer

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = REPO_ROOT / "analysis" / "la_snp_1kg_frequencies.csv"
DEFAULT_OUTPUT = REPO_ROOT / "analysis" / "la_snp_af_1kg_vs_gnomad.csv"
DEFAULT_SCATTER = REPO_ROOT / "figures" / "af_1kg_vs_gnomad_scatter.png"
DEFAULT_SUMMARY = REPO_ROOT / "analysis" / "af_comparison_summary.md"
DEFAULT_CACHE = REPO_ROOT / "data" / "geo" / "gnomad_r4_nfe_cache.json"

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET = "gnomad_r4"
GNOMAD_POPULATION = "nfe"
REFERENCE_GENOME = "GRCh38"

DEFAULT_BATCH_SIZE = 8
DEFAULT_MIN_INTERVAL_SEC = 0.75
DEFAULT_TIMEOUT_SEC = 45.0
DEFAULT_MAX_RETRIES = 4
DIFF_THRESHOLD = 0.05

LOG = logging.getLogger(__name__)

LogLevel = Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
_LOG_LEVEL_CHOICE = click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

VARIANT_FIELDS = """
  variant_id
  rsids
  chrom
  pos
  ref
  alt
  joint { populations { id ac an } }
  exome { populations { id ac an } }
  genome { populations { id ac an } }
"""

app = typer.Typer(
    add_completion=False,
    help="Fetch gnomAD v4 NFE AFs and compare to 1KG frequencies, or summarize a comparison CSV.",
)


def normalize_chromosome_label(chromosome: str) -> str:
    """Normalise a chromosome label for gnomAD variant IDs (strip ``chr``, map MT).

    Args:
        chromosome: Raw chromosome label.

    Returns:
        Compact label without a ``chr`` prefix (``MT`` for mitochondrial).
    """
    label = str(chromosome).strip()
    if label.upper().startswith("CHR"):
        label = label[3:]
    numeric = pd.to_numeric(label, errors="coerce")
    if pd.notna(numeric) and float(numeric).is_integer():
        label = str(int(numeric))
    upper = label.upper()
    if upper in {"M", "MT"}:
        return "MT"
    return upper if upper in {"X", "Y", "MT"} else label


def normalize_rsid(rsid: str) -> str:
    """Ensure an rsID string has an ``rs`` prefix.

    Args:
        rsid: Candidate rs identifier (with or without ``rs``).

    Returns:
        Normalised rsID string, or ``""`` if ``rsid`` is blank.
    """
    value = str(rsid).strip()
    if not value:
        return ""
    return value if value.lower().startswith("rs") else f"rs{value}"


def to_gnomad_variant_id(chromosome: str, position: int, ref: str, alt: str) -> str:
    """Build a gnomAD-style variant ID ``CHROM-POS-REF-ALT``.

    Args:
        chromosome: Chromosome label.
        position: One-based genomic position.
        ref: Reference allele.
        alt: Alternate allele.

    Returns:
        Variant ID string suitable for the gnomAD GraphQL ``variantId`` field.
    """
    return f"{normalize_chromosome_label(chromosome)}-{position}-{ref}-{alt}"


def read_1kg_frequencies(path: Path) -> pd.DataFrame:
    """Load the 1KG allele-frequency table produced by ``rogen-ukb-manifest extract``.

    Args:
        path: Path to the 1KG frequency CSV.

    Returns:
        DataFrame including at least ``rsID`` and ``AF``.

    Raises:
        ValueError: If required columns are missing.
    """
    df = pd.read_csv(path)
    required = {"rsID", "AF"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {sorted(missing)}")
    return df


def load_cache(path: Path) -> dict[str, dict[str, Any]]:
    """Load a gnomAD GraphQL response cache from JSON, or return an empty dict.

    Args:
        path: JSON cache file path.

    Returns:
        Mapping from rsID to cached lookup payloads.

    Raises:
        ValueError: If the file exists but is not a JSON object.
    """
    if not path.is_file():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Cache file must contain a JSON object: {path}")
    return payload


def save_cache(path: Path, cache: dict[str, dict[str, Any]]) -> None:
    """Persist the gnomAD GraphQL response cache as formatted JSON.

    Args:
        path: Destination JSON path; parent directories are created.
        cache: Mapping from rsID to cached lookup payloads.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cache, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def population_nfe_af(variant: dict[str, Any] | None) -> float | None:
    """Extract the non-Finnish European allele frequency from a gnomAD variant payload.

    Args:
        variant: gnomAD GraphQL variant object, or ``None``.

    Returns:
        NFE allele frequency from joint/exome/genome populations, or ``None``
        if unavailable.
    """
    if not variant:
        return None
    for source in ("joint", "exome", "genome"):
        block = variant.get(source)
        if not isinstance(block, dict):
            continue
        populations = block.get("populations")
        if not isinstance(populations, list):
            continue
        for population in populations:
            if population.get("id") != GNOMAD_POPULATION:
                continue
            allele_number = population.get("an")
            allele_count = population.get("ac")
            if allele_number and allele_number > 0 and allele_count is not None:
                return float(allele_count) / float(allele_number)
    return None


def cache_entry_from_variant(
    variant: dict[str, Any] | None,
    *,
    lookup_method: str,
) -> dict[str, Any]:
    """Build a cache record from a gnomAD variant payload.

    Args:
        variant: gnomAD GraphQL variant object, or ``None`` when unresolved.
        lookup_method: Label for how the lookup was performed (e.g.
            ``variant_id``, ``region_rsid``, ``unresolved``).

    Returns:
        Cache entry with ``af_gnomad_nfe``, ``variant_id``, ``lookup_method``,
        and ``fetched_at``.
    """
    return {
        "af_gnomad_nfe": population_nfe_af(variant),
        "variant_id": None if not variant else variant.get("variant_id"),
        "lookup_method": lookup_method,
        "fetched_at": datetime.now(tz=UTC).isoformat(),
    }


class GnomadClient:
    """Polite gnomAD GraphQL client with retry and pacing."""

    def __init__(
        self,
        *,
        timeout_sec: float,
        min_interval_sec: float,
        max_retries: int,
        session: requests.Session | None = None,
    ) -> None:
        """Initialise the GraphQL client.

        Args:
            timeout_sec: Per-request socket timeout.
            min_interval_sec: Minimum seconds between completed requests.
            max_retries: Maximum retries for transient HTTP/network failures.
            session: Optional :class:`requests.Session` for connection pooling.
        """
        self.timeout_sec = timeout_sec
        self.min_interval_sec = min_interval_sec
        self.max_retries = max_retries
        self._session = session or requests.Session()
        self._last_request_end: float | None = None
        self._session.headers.setdefault(
            "User-Agent",
            "rogen-aging-compare-af-gnomad/0.1 (ROGEN; academic research; gnomAD GraphQL)",
        )
        self._session.headers.setdefault("Content-Type", "application/json")

    def close(self) -> None:
        """Close the underlying HTTP session."""
        self._session.close()

    def _pace(self) -> None:
        if self._last_request_end is None:
            return
        elapsed = time.monotonic() - self._last_request_end
        wait = self.min_interval_sec - elapsed
        if wait > 0:
            time.sleep(wait)

    def _post(self, query: str, variables: dict[str, Any] | None = None) -> dict[str, Any]:
        self._pace()
        payload: dict[str, Any] = {"query": query}
        if variables is not None:
            payload["variables"] = variables

        attempt = 0
        while True:
            attempt += 1
            try:
                response = self._session.post(GNOMAD_API_URL, json=payload, timeout=self.timeout_sec)
            except requests.RequestException as exc:
                if attempt > self.max_retries:
                    raise RuntimeError(f"gnomAD request failed after {self.max_retries} attempts: {exc}") from exc
                sleep_s = self.min_interval_sec * (2 ** (attempt - 1))
                LOG.warning("gnomAD request error (%s); retrying in %.2fs", exc, sleep_s)
                time.sleep(sleep_s)
                continue

            self._last_request_end = time.monotonic()

            if response.status_code in {429, 502, 503, 504}:
                if attempt > self.max_retries:
                    response.raise_for_status()
                retry_after = response.headers.get("Retry-After")
                sleep_s = float(retry_after) if retry_after else self.min_interval_sec * (2 ** (attempt - 1))
                LOG.warning("gnomAD HTTP %s; sleeping %.2fs", response.status_code, sleep_s)
                time.sleep(sleep_s)
                continue

            response.raise_for_status()
            body = response.json()
            if "errors" in body and body.get("data") in (None, {}):
                messages = "; ".join(error.get("message", "unknown error") for error in body["errors"])
                raise RuntimeError(f"gnomAD GraphQL error: {messages}")
            return body

    def fetch_variants_by_id(self, variant_ids: list[str]) -> dict[str, dict[str, Any] | None]:
        """Fetch gnomAD variants by ``chrom-pos-ref-alt`` IDs in one GraphQL query.

        Args:
            variant_ids: gnomAD variant ID strings.

        Returns:
            Mapping from requested variant ID to the GraphQL variant payload
            (or ``None`` when not found).

        Raises:
            RuntimeError: On exhausted retries or GraphQL-level failures.
        """
        if not variant_ids:
            return {}

        alias_blocks: list[str] = []
        for index, variant_id in enumerate(variant_ids):
            safe_id = variant_id.replace('"', "")
            alias_blocks.append(
                f'v{index}: variant(variantId: "{safe_id}", dataset: {GNOMAD_DATASET}) {{ {VARIANT_FIELDS} }}'
            )
        query = "query {\n" + "\n".join(alias_blocks) + "\n}"
        body = self._post(query)
        data = body.get("data") or {}

        results: dict[str, dict[str, Any] | None] = {}
        for index, variant_id in enumerate(variant_ids):
            results[variant_id] = data.get(f"v{index}")
        return results

    def fetch_variant_by_region(
        self,
        chromosome: str,
        position: int,
        rsid: str,
    ) -> dict[str, Any] | None:
        """Look up a gnomAD variant by single-base region and match on rsID.

        Args:
            chromosome: Chromosome label.
            position: One-based genomic position.
            rsid: Expected rs identifier among region hits.

        Returns:
            Matching GraphQL variant payload, or ``None`` if not found.

        Raises:
            RuntimeError: On exhausted retries or GraphQL-level failures.
        """
        chrom = normalize_chromosome_label(chromosome)
        query = """
        query RegionLookup($chrom: String!, $start: Int!, $stop: Int!, $referenceGenome: ReferenceGenomeId!) {
          region(chrom: $chrom, start: $start, stop: $stop, reference_genome: $referenceGenome) {
            variants(dataset: gnomad_r4) {
              variant_id
              rsids
              ref
              alt
              joint { populations { id ac an } }
              exome { populations { id ac an } }
              genome { populations { id ac an } }
            }
          }
        }
        """
        variables = {
            "chrom": chrom,
            "start": position,
            "stop": position,
            "referenceGenome": REFERENCE_GENOME,
        }
        body = self._post(query, variables)
        region = (body.get("data") or {}).get("region")
        if not region:
            return None
        variants = region.get("variants") or []
        target = normalize_rsid(rsid).lower()
        for variant in variants:
            rsids = variant.get("rsids") or []
            normalized = {normalize_rsid(value).lower() for value in rsids}
            if target in normalized:
                return variant
        return None


@dataclass
class LookupPlan:
    """One LA-SNP row scheduled for a gnomAD allele-frequency lookup."""

    rsid: str
    af_1kg: float | None
    variant_id: str | None
    chromosome: str | None
    position: int | None
    ref: str | None
    alt: str | None


def build_lookup_plans(df: pd.DataFrame) -> list[LookupPlan]:
    """Build gnomAD lookup plans from a 1KG frequency table with optional coordinates.

    Args:
        df: 1KG frequency table with ``rsID``, ``AF``, and optionally
            ``chrom``, ``pos``, ``ref``, ``alt``.

    Returns:
        One :class:`LookupPlan` per non-empty rsID row.
    """
    plans: list[LookupPlan] = []
    for _, row in df.iterrows():
        rsid = normalize_rsid(row["rsID"]) if pd.notna(row["rsID"]) else ""
        if not rsid:
            continue

        af_1kg = None if pd.isna(row["AF"]) else float(row["AF"])
        chrom = None if pd.isna(row.get("chrom")) else normalize_chromosome_label(str(row["chrom"]))
        position = None
        if pd.notna(row.get("pos")):
            position = int(float(row["pos"]))

        ref = "" if pd.isna(row.get("ref")) else str(row["ref"]).strip()
        alt = "" if pd.isna(row.get("alt")) else str(row["alt"]).strip()

        variant_id = None
        if chrom and position and ref and alt:
            variant_id = to_gnomad_variant_id(chrom, position, ref, alt)

        plans.append(
            LookupPlan(
                rsid=rsid,
                af_1kg=af_1kg,
                variant_id=variant_id,
                chromosome=chrom,
                position=position,
                ref=ref or None,
                alt=alt or None,
            )
        )
    return plans


def fetch_gnomad_afs(
    plans: list[LookupPlan],
    *,
    cache_path: Path,
    refresh_cache: bool,
    batch_size: int,
    client: GnomadClient,
) -> dict[str, dict[str, Any]]:
    """Fetch gnomAD allele frequencies for lookup plans, using and updating ``cache_path``.

    Args:
        plans: Per-SNP lookup plans from :func:`build_lookup_plans`.
        cache_path: JSON cache path for rsID → AF payloads.
        refresh_cache: If true, ignore existing cache entries and re-query.
        batch_size: Number of variant IDs per GraphQL request.
        client: Configured :class:`GnomadClient`.

    Returns:
        Updated cache mapping from rsID to lookup payloads.
    """
    cache = {} if refresh_cache else load_cache(cache_path)
    pending_variant_ids: dict[str, list[LookupPlan]] = {}
    pending_regions: dict[str, LookupPlan] = {}

    for plan in plans:
        cached = cache.get(plan.rsid)
        if cached is not None and not refresh_cache:
            continue
        if plan.variant_id:
            pending_variant_ids.setdefault(plan.variant_id, []).append(plan)
        elif plan.chromosome and plan.position is not None:
            pending_regions[plan.rsid] = plan
        else:
            cache[plan.rsid] = cache_entry_from_variant(None, lookup_method="unresolved")
            LOG.warning("Cannot resolve gnomAD lookup for %s (missing coordinates).", plan.rsid)

    variant_id_list = list(pending_variant_ids)
    for batch_start in range(0, len(variant_id_list), batch_size):
        batch_ids = variant_id_list[batch_start : batch_start + batch_size]
        LOG.info(
            "Querying gnomAD for variant batch %s-%s of %s",
            batch_start + 1,
            batch_start + len(batch_ids),
            len(variant_id_list),
        )
        batch_results = client.fetch_variants_by_id(batch_ids)
        for variant_id in batch_ids:
            variant = batch_results.get(variant_id)
            for plan in pending_variant_ids[variant_id]:
                cache[plan.rsid] = cache_entry_from_variant(variant, lookup_method="variant_id")
                if variant is None:
                    pending_regions[plan.rsid] = plan
                elif cache[plan.rsid]["af_gnomad_nfe"] is None:
                    LOG.warning("SNP %s found in gnomAD but lacks NFE AF.", plan.rsid)

    for plan in pending_regions.values():
        if cache.get(plan.rsid, {}).get("af_gnomad_nfe") is not None:
            continue
        if not plan.chromosome or plan.position is None:
            continue
        LOG.info("Region fallback lookup for %s at %s:%s", plan.rsid, plan.chromosome, plan.position)
        variant = client.fetch_variant_by_region(plan.chromosome, plan.position, plan.rsid)
        cache[plan.rsid] = cache_entry_from_variant(variant, lookup_method="region_rsid")
        if variant is None:
            LOG.warning("SNP %s not found in gnomAD v4 (%s).", plan.rsid, GNOMAD_DATASET)

    save_cache(cache_path, cache)
    return cache


def build_comparison_table(plans: list[LookupPlan], cache: dict[str, dict[str, Any]]) -> pd.DataFrame:
    """Merge lookup plans with cached gnomAD responses into a comparison table.

    Args:
        plans: Per-SNP lookup plans.
        cache: rsID → cached gnomAD lookup payloads.

    Returns:
        Comparison DataFrame with ``rsID``, ``AF_1kg``, ``AF_gnomad_nfe``,
        ``abs_diff``, and ``large_diff``.
    """
    rows: list[dict[str, Any]] = []
    for plan in plans:
        cached = cache.get(plan.rsid, {})
        af_gnomad = cached.get("af_gnomad_nfe")
        af_1kg = plan.af_1kg
        abs_diff = None
        large_diff = False
        if af_1kg is not None and af_gnomad is not None:
            abs_diff = abs(af_1kg - float(af_gnomad))
            large_diff = abs_diff > DIFF_THRESHOLD

        rows.append(
            {
                "rsID": plan.rsid,
                "AF_1kg": af_1kg,
                "AF_gnomad_nfe": af_gnomad,
                "abs_diff": abs_diff,
                "large_diff": large_diff,
            }
        )
    return pd.DataFrame(rows)


def plot_scatter(comparison: pd.DataFrame, output_path: Path) -> None:
    """Write a 1KG vs gnomAD NFE allele-frequency scatter plot.

    Args:
        comparison: Output of :func:`build_comparison_table`.
        output_path: Destination PNG path.
    """
    plotted = comparison.dropna(subset=["AF_1kg", "AF_gnomad_nfe"])
    if plotted.empty:
        LOG.warning("No overlapping AF values to plot; skipping scatter.")
        return

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    normal = plotted[~plotted["large_diff"]]
    flagged = plotted[plotted["large_diff"]]

    ax.scatter(normal["AF_1kg"], normal["AF_gnomad_nfe"], s=36, color="#5B7C99", alpha=0.85, label="|diff| ≤ 0.05")
    if not flagged.empty:
        ax.scatter(
            flagged["AF_1kg"],
            flagged["AF_gnomad_nfe"],
            s=48,
            color="#C45C3E",
            alpha=0.95,
            label="|diff| > 0.05",
        )

    max_af = max(plotted["AF_1kg"].max(), plotted["AF_gnomad_nfe"].max())
    pad = 0.02
    upper = min(1.0, max_af + pad)
    ax.plot([0, upper], [0, upper], linestyle="--", color="#666666", linewidth=1, label="Identity")
    ax.set_xlim(0, upper)
    ax.set_ylim(0, upper)
    ax.set_xlabel("AF (1000 Genomes)")
    ax.set_ylabel("AF (gnomAD v4, NFE)")
    ax.set_title("LA-SNP allele frequencies: 1KG vs gnomAD NFE")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="lower right", frameon=False, fontsize=8)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    LOG.info("Wrote scatter plot to %s", output_path.resolve())


@dataclass(frozen=True)
class AfComparisonSummary:
    """Headline statistics from a 1KG vs gnomAD allele-frequency comparison."""

    total_snps: int
    resolved_both: int
    missing_1kg: int
    missing_gnomad: int
    concordant: int
    discordant: int
    mean_abs_diff_concordant: float | None
    median_abs_diff_concordant: float | None
    top_discordant: pd.DataFrame
    diff_threshold: float


def read_comparison_table(path: Path) -> pd.DataFrame:
    """Load a comparison CSV written by ``compare_main``.

    Args:
        path: Path to the comparison CSV.

    Returns:
        DataFrame with at least ``rsID``, ``AF_1kg``, ``AF_gnomad_nfe``, and
        ``abs_diff`` (computed if missing).

    Raises:
        ValueError: If required columns are missing.
    """
    df = pd.read_csv(path)
    required = {"rsID", "AF_1kg", "AF_gnomad_nfe"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Comparison CSV is missing required columns: {sorted(missing)}")
    if "abs_diff" not in df.columns:
        paired = df["AF_1kg"].notna() & df["AF_gnomad_nfe"].notna()
        df = df.copy()
        df["abs_diff"] = None
        df.loc[paired, "abs_diff"] = (df.loc[paired, "AF_1kg"] - df.loc[paired, "AF_gnomad_nfe"]).abs()
    return df


def summarize_comparison(
    comparison: pd.DataFrame,
    *,
    diff_threshold: float = DIFF_THRESHOLD,
) -> AfComparisonSummary:
    """Compute headline concordance statistics for a 1KG vs gnomAD table.

    Args:
        comparison: Comparison table from :func:`build_comparison_table` or
            :func:`read_comparison_table`.
        diff_threshold: Absolute AF difference cutoff; loci with
            ``|ΔAF| < threshold`` are treated as concordant.

    Returns:
        :class:`AfComparisonSummary` with counts and top discordant loci.
    """
    total_snps = len(comparison)
    has_1kg = comparison["AF_1kg"].notna()
    has_gnomad = comparison["AF_gnomad_nfe"].notna()
    resolved_both = int((has_1kg & has_gnomad).sum())
    missing_1kg = int((~has_1kg).sum())
    missing_gnomad = int((~has_gnomad).sum())

    paired = comparison.loc[has_1kg & has_gnomad].copy()
    paired["abs_diff"] = (paired["AF_1kg"] - paired["AF_gnomad_nfe"]).abs()
    concordant_mask = paired["abs_diff"] < diff_threshold
    concordant = int(concordant_mask.sum())
    discordant = int((~concordant_mask).sum())

    concordant_diffs = paired.loc[concordant_mask, "abs_diff"]
    mean_abs_diff = float(concordant_diffs.mean()) if not concordant_diffs.empty else None
    median_abs_diff = float(concordant_diffs.median()) if not concordant_diffs.empty else None

    top_discordant = (
        paired.sort_values("abs_diff", ascending=False)
        .head(5)[["rsID", "AF_1kg", "AF_gnomad_nfe", "abs_diff"]]
        .reset_index(drop=True)
    )

    return AfComparisonSummary(
        total_snps=total_snps,
        resolved_both=resolved_both,
        missing_1kg=missing_1kg,
        missing_gnomad=missing_gnomad,
        concordant=concordant,
        discordant=discordant,
        mean_abs_diff_concordant=mean_abs_diff,
        median_abs_diff_concordant=median_abs_diff,
        top_discordant=top_discordant,
        diff_threshold=diff_threshold,
    )


def _format_af(value: float | None) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):.4f}"


def format_summary_markdown(summary: AfComparisonSummary) -> str:
    """Render an :class:`AfComparisonSummary` as Markdown prose plus a table.

    Args:
        summary: Output of :func:`summarize_comparison`.

    Returns:
        Markdown string with a summary paragraph and top-discordant table.
    """
    threshold = summary.diff_threshold
    mean_part = (
        "N/A"
        if summary.mean_abs_diff_concordant is None
        else f"{summary.mean_abs_diff_concordant:.4f}"
    )
    median_part = (
        "N/A"
        if summary.median_abs_diff_concordant is None
        else f"{summary.median_abs_diff_concordant:.4f}"
    )

    paragraph = (
        f"Of {summary.total_snps} LA-SNPs in the comparison table, "
        f"{summary.resolved_both} had allele frequencies in both 1000 Genomes and gnomAD v4 NFE; "
        f"{summary.missing_1kg} were missing from 1KG and {summary.missing_gnomad} from gnomAD. "
        f"Among the {summary.resolved_both} with paired frequencies, "
        f"{summary.concordant} were concordant (|ΔAF| < {threshold:g}) and "
        f"{summary.discordant} were discordant (|ΔAF| ≥ {threshold:g}). "
        f"Mean and median |ΔAF| among concordant loci were {mean_part} and {median_part}, respectively."
    )

    table_lines = [
        "",
        "| rsID | AF_1kg | AF_gnomad_nfe | abs_diff |",
        "| --- | ---: | ---: | ---: |",
    ]
    if summary.top_discordant.empty:
        table_lines.append("| — | — | — | — |")
    else:
        for _, row in summary.top_discordant.iterrows():
            table_lines.append(
                "| {rsid} | {af_1kg} | {af_gnomad} | {abs_diff} |".format(
                    rsid=row["rsID"],
                    af_1kg=_format_af(row["AF_1kg"]),
                    af_gnomad=_format_af(row["AF_gnomad_nfe"]),
                    abs_diff=_format_af(row["abs_diff"]),
                )
            )

    return paragraph + "\n".join(table_lines) + "\n"


def log_missing_gnomad(comparison: pd.DataFrame) -> None:
    """Log rsIDs missing gnomAD NFE allele frequencies.

    Args:
        comparison: Comparison table with ``rsID`` and ``AF_gnomad_nfe``.
    """
    missing = comparison.loc[comparison["AF_gnomad_nfe"].isna(), "rsID"].tolist()
    if not missing:
        LOG.info("All %s SNPs were found in gnomAD v4 with NFE frequencies.", len(comparison))
        return
    LOG.warning(
        "%s SNP(s) missing from gnomAD v4 or lacking NFE AF: %s",
        len(missing),
        ", ".join(missing),
    )


def compare_main(
    *,
    input_path: Path = DEFAULT_INPUT,
    output: Path = DEFAULT_OUTPUT,
    scatter: Path = DEFAULT_SCATTER,
    cache: Path = DEFAULT_CACHE,
    refresh_cache: bool = False,
    batch_size: int = DEFAULT_BATCH_SIZE,
    min_interval: float = DEFAULT_MIN_INTERVAL_SEC,
    timeout: float = DEFAULT_TIMEOUT_SEC,
    max_retries: int = DEFAULT_MAX_RETRIES,
    log_level: LogLevel = "INFO",
) -> int:
    """Fetch gnomAD AFs and write comparison CSV plus scatter plot.

    Args:
        input_path: 1KG frequency CSV path.
        output: Destination comparison CSV path.
        scatter: Destination scatter PNG path.
        cache: JSON cache path for gnomAD responses.
        refresh_cache: If true, ignore cached rsIDs and re-query.
        batch_size: Variants per GraphQL request.
        min_interval: Minimum seconds between gnomAD HTTP requests.
        timeout: HTTP timeout in seconds.
        max_retries: Retries per HTTP request.
        log_level: Logging verbosity name.

    Returns:
        Process exit code (``0`` on success, ``1`` on failure).
    """
    logging.basicConfig(
        level=getattr(logging, str(log_level)),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    if not input_path.is_file():
        LOG.error("Input CSV does not exist: %s", input_path.resolve())
        return 1

    df = read_1kg_frequencies(input_path)
    plans = build_lookup_plans(df)
    LOG.info("Loaded %s SNP rows from %s", len(plans), input_path.resolve())

    client = GnomadClient(
        timeout_sec=float(timeout),
        min_interval_sec=float(min_interval),
        max_retries=int(max_retries),
    )
    try:
        cache_data = fetch_gnomad_afs(
            plans,
            cache_path=cache,
            refresh_cache=bool(refresh_cache),
            batch_size=int(batch_size),
            client=client,
        )
    finally:
        client.close()

    comparison = build_comparison_table(plans, cache_data)
    log_missing_gnomad(comparison)

    output.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(output, index=False)
    LOG.info("Wrote comparison table (%s rows) to %s", len(comparison), output.resolve())

    plot_scatter(comparison, scatter)
    return 0


def summarize_main(
    *,
    input_path: Path = DEFAULT_OUTPUT,
    output: Path = DEFAULT_SUMMARY,
    diff_threshold: float = DIFF_THRESHOLD,
    log_level: LogLevel = "INFO",
) -> int:
    """Summarize an existing comparison CSV as Markdown.

    Args:
        input_path: Comparison CSV from the compare step.
        output: Destination Markdown summary path.
        diff_threshold: Concordance cutoff for ``|ΔAF|``.
        log_level: Logging verbosity name.

    Returns:
        Process exit code (``0`` on success, ``1`` on failure).
    """
    logging.basicConfig(
        level=getattr(logging, str(log_level)),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    if not input_path.is_file():
        LOG.error("Comparison CSV does not exist: %s", input_path.resolve())
        return 1

    comparison = read_comparison_table(input_path)
    summary = summarize_comparison(comparison, diff_threshold=float(diff_threshold))
    markdown = format_summary_markdown(summary)

    print(markdown, end="")

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(markdown, encoding="utf-8")
    LOG.info("Wrote summary to %s", output.resolve())
    return 0


@app.command("compare")
def compare_cmd(
    input_path: Path = typer.Option(DEFAULT_INPUT, "--input", path_type=Path, help="1KG frequency CSV"),
    output: Path = typer.Option(DEFAULT_OUTPUT, "--output", path_type=Path, help="Comparison CSV"),
    scatter: Path = typer.Option(DEFAULT_SCATTER, "--scatter", path_type=Path, help="Scatter PNG path"),
    cache: Path = typer.Option(DEFAULT_CACHE, "--cache", path_type=Path, help="JSON cache path"),
    refresh_cache: bool = typer.Option(False, "--refresh-cache", help="Ignore cached rsIDs and re-query gnomAD"),
    batch_size: int = typer.Option(DEFAULT_BATCH_SIZE, "--batch-size", help="Variants per GraphQL request"),
    min_interval: float = typer.Option(
        DEFAULT_MIN_INTERVAL_SEC,
        "--min-interval",
        help="Minimum seconds between gnomAD HTTP requests",
    ),
    timeout: float = typer.Option(DEFAULT_TIMEOUT_SEC, "--timeout", help="HTTP timeout in seconds"),
    max_retries: int = typer.Option(DEFAULT_MAX_RETRIES, "--max-retries", help="Retries per HTTP request"),
    log_level: str = typer.Option(
        "INFO",
        "--log-level",
        click_type=_LOG_LEVEL_CHOICE,
        help="Logging verbosity",
    ),
) -> int:
    """Fetch gnomAD v4 NFE AFs and compare to 1KG frequencies."""
    return compare_main(
        input_path=input_path,
        output=output,
        scatter=scatter,
        cache=cache,
        refresh_cache=refresh_cache,
        batch_size=batch_size,
        min_interval=min_interval,
        timeout=timeout,
        max_retries=max_retries,
        log_level=cast(LogLevel, log_level),
    )


@app.command("summarize")
def summarize_cmd(
    input_path: Path = typer.Option(
        DEFAULT_OUTPUT,
        "--input",
        path_type=Path,
        help="Comparison CSV from the compare step",
    ),
    output: Path = typer.Option(
        DEFAULT_SUMMARY,
        "--output",
        path_type=Path,
        help="Markdown summary path",
    ),
    diff_threshold: float = typer.Option(
        DIFF_THRESHOLD,
        "--diff-threshold",
        help="Concordance cutoff; discordant when |ΔAF| is at or above this value",
    ),
    log_level: str = typer.Option(
        "INFO",
        "--log-level",
        click_type=_LOG_LEVEL_CHOICE,
        help="Logging verbosity",
    ),
) -> int:
    """Summarize a 1KG vs gnomAD comparison CSV for reporting."""
    return summarize_main(
        input_path=input_path,
        output=output,
        diff_threshold=diff_threshold,
        log_level=cast(LogLevel, log_level),
    )


def main(argv: list[str] | None = None) -> int:
    """Dispatch to compare (default) or summarize subcommands.

    Args:
        argv: Argument vector without the program name. Defaults to
            ``sys.argv[1:]``. Bare flags default to the ``compare`` command.

    Returns:
        Process exit code from the selected Typer command.
    """
    if argv is None:
        argv = sys.argv[1:]
    if not argv or argv[0].startswith("-"):
        dispatch = ["compare", *argv]
    elif argv[0] in {"compare", "summarize"}:
        dispatch = list(argv)
    else:
        dispatch = ["compare", *argv]

    try:
        result = app(args=dispatch, standalone_mode=False)
    except typer.Exit as exc:
        return int(exc.exit_code)
    except SystemExit as exc:
        code = exc.code
        if code is None:
            return 0
        if isinstance(code, int):
            return code
        return 1
    return 0 if result is None else int(result)


if __name__ == "__main__":
    sys.exit(main())
