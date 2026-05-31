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
        --scatter analysis/af_1kg_vs_gnomad_scatter.png
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
import requests

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = REPO_ROOT / "analysis" / "la_snp_1kg_frequencies.csv"
DEFAULT_OUTPUT = REPO_ROOT / "analysis" / "la_snp_af_1kg_vs_gnomad.csv"
DEFAULT_SCATTER = REPO_ROOT / "analysis" / "af_1kg_vs_gnomad_scatter.png"
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


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="1KG frequency CSV")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Comparison CSV")
    parser.add_argument("--scatter", type=Path, default=DEFAULT_SCATTER, help="Scatter PNG path")
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE, help="JSON cache path")
    parser.add_argument(
        "--refresh-cache",
        action="store_true",
        help="Ignore cached rsIDs and re-query gnomAD",
    )
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE, help="Variants per GraphQL request")
    parser.add_argument(
        "--min-interval",
        type=float,
        default=DEFAULT_MIN_INTERVAL_SEC,
        help="Minimum seconds between gnomAD HTTP requests",
    )
    parser.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT_SEC, help="HTTP timeout in seconds")
    parser.add_argument("--max-retries", type=int, default=DEFAULT_MAX_RETRIES, help="Retries per HTTP request")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging verbosity",
    )
    return parser.parse_args(argv)


def normalize_chromosome_label(chromosome: str) -> str:
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
    value = str(rsid).strip()
    if not value:
        return ""
    return value if value.lower().startswith("rs") else f"rs{value}"


def to_gnomad_variant_id(chromosome: str, position: int, ref: str, alt: str) -> str:
    return f"{normalize_chromosome_label(chromosome)}-{position}-{ref}-{alt}"


def read_1kg_frequencies(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"rsID", "AF"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {sorted(missing)}")
    return df


def load_cache(path: Path) -> dict[str, dict[str, Any]]:
    if not path.is_file():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Cache file must contain a JSON object: {path}")
    return payload


def save_cache(path: Path, cache: dict[str, dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cache, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def population_nfe_af(variant: dict[str, Any] | None) -> float | None:
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
    rsid: str
    af_1kg: float | None
    variant_id: str | None
    chromosome: str | None
    position: int | None
    ref: str | None
    alt: str | None


def build_lookup_plans(df: pd.DataFrame) -> list[LookupPlan]:
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


def log_missing_gnomad(comparison: pd.DataFrame) -> None:
    missing = comparison.loc[comparison["AF_gnomad_nfe"].isna(), "rsID"].tolist()
    if not missing:
        LOG.info("All %s SNPs were found in gnomAD v4 with NFE frequencies.", len(comparison))
        return
    LOG.warning(
        "%s SNP(s) missing from gnomAD v4 or lacking NFE AF: %s",
        len(missing),
        ", ".join(missing),
    )


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, str(args.log_level)),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    input_path: Path = args.input
    if not input_path.is_file():
        LOG.error("Input CSV does not exist: %s", input_path.resolve())
        return 1

    df = read_1kg_frequencies(input_path)
    plans = build_lookup_plans(df)
    LOG.info("Loaded %s SNP rows from %s", len(plans), input_path.resolve())

    client = GnomadClient(
        timeout_sec=float(args.timeout),
        min_interval_sec=float(args.min_interval),
        max_retries=int(args.max_retries),
    )
    try:
        cache = fetch_gnomad_afs(
            plans,
            cache_path=args.cache,
            refresh_cache=bool(args.refresh_cache),
            batch_size=int(args.batch_size),
            client=client,
        )
    finally:
        client.close()

    comparison = build_comparison_table(plans, cache)
    log_missing_gnomad(comparison)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(args.output, index=False)
    LOG.info("Wrote comparison table (%s rows) to %s", len(comparison), args.output.resolve())

    plot_scatter(comparison, args.scatter)
    return 0


if __name__ == "__main__":
    sys.exit(main())
