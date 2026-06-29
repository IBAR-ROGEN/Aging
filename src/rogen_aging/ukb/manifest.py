"""Offline UK Biobank (UKB) LA-SNP manifest builder and public-frequency proxy for ROGEN.

v0.1 — **build** subcommand (default when no subcommand is given): reads a
gene–SNP overlap table from Excel, resolves each dbSNP rs identifier to a
GRCh38 chromosomal locus via the Ensembl Variation REST API, and writes a CSV
manifest suitable for downstream bulk genotype extraction (for example,
imputed genotype resources often discussed in conjunction with UK Biobank
showcase field **22418** — genotype/imputed bulk files).

v0.2 — **extract** subcommand: uses that manifest to pull the same SNP loci
from publicly available 1000 Genomes Project GRCh38 VCFs (bgzipped and
tabix-indexed) and reports per-variant allele frequency and call count. This
is a **public-data proxy** for expected UKB allele frequencies during
development; it does **not** access UK Biobank participant genotypes and
makes **no** DNAnexus, dx-toolkit, or dxFUSE calls.

The manifest encodes chromosome and GRCh38 position for harmonisation with
contemporary reference-based workflows. Native UKB imputed BGEN products
are historically GRCh37-aligned; confirm liftover requirements before
matching coordinates to released UKB variant sets.

Examples:
    uv run rogen-ukb-manifest build \\
        --input overlapping_genes_with_snps.xlsx \\
        --output analysis/ukb_snp_manifest_v0.1.csv

    uv run rogen-ukb-manifest extract \\
        --manifest analysis/ukb_snp_manifest_v0.1.csv \\
        --vcf-glob 'data/1kg/ALL.chr*.vcf.gz' \\
        --output analysis/la_snp_1kg_frequencies.csv
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import re
import sys
import time
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import quote

import cyvcf2
import pandas as pd
import requests

ENSEMBL_VARIATION_BASE = "https://rest.ensembl.org/variation/human"
DEFAULT_ASSEMBLY = "GRCh38"
# Conservative pacing; Ensembl recommends modest parallelism and polite usage.
DEFAULT_MIN_INTERVAL_SEC = 0.34
DEFAULT_TIMEOUT_SEC = 30.0
DEFAULT_MAX_RETRIES = 4
DEFAULT_MANIFEST_CSV = Path("analysis/ukb_snp_manifest_v0.1.csv")
DEFAULT_1KG_FREQ_CSV = Path("analysis/la_snp_1kg_frequencies.csv")
_VCF_GLOB_SUFFIXES = ("*.vcf.gz", "*.vcf.bgz", "*.bcf")
_CHROM_IN_FILENAME_RE = re.compile(
    r"(?:^|[_.-])((?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT|M))(?:[_.-]|\.|$)",
    re.IGNORECASE,
)

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
        "rogen-aging-ukb-manifest/0.2 (ROGEN; academic research; Ensembl REST)",
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


def normalize_chromosome_label(chromosome: str) -> str:
    """Normalise a chromosome label to a compact form (no ``chr`` prefix)."""
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


def chromosome_query_aliases(chromosome: str) -> list[str]:
    """Return plausible contig names for tabix region queries."""
    base = normalize_chromosome_label(chromosome)
    aliases = [base, f"chr{base}"]
    if base == "MT":
        aliases.extend(["M", "chrM"])
    seen: set[str] = set()
    ordered: list[str] = []
    for alias in aliases:
        if alias not in seen:
            seen.add(alias)
            ordered.append(alias)
    return ordered


def expand_vcf_paths(vcf_glob: str) -> list[Path]:
    """Expand a filesystem path or glob to indexed 1KG-style VCF paths."""
    matches = sorted(Path(p) for p in glob.glob(vcf_glob))
    if matches:
        return matches

    candidate = Path(vcf_glob)
    if candidate.is_dir():
        paths: list[Path] = []
        for pattern in _VCF_GLOB_SUFFIXES:
            paths.extend(sorted(candidate.glob(pattern)))
        return paths

    if candidate.is_file():
        return [candidate]

    raise FileNotFoundError(f"No VCF files matched: {vcf_glob}")


def infer_vcf_chromosomes(path: Path) -> set[str]:
    """Infer chromosome labels encoded in a VCF filename, if any."""
    labels: set[str] = set()
    for match in _CHROM_IN_FILENAME_RE.findall(path.name):
        labels.add(normalize_chromosome_label(match))
    return labels


def build_chromosome_vcf_index(vcf_paths: Sequence[Path]) -> tuple[dict[str, Path], Path | None]:
    """Map normalised chromosome labels to per-chromosome VCF paths."""
    index: dict[str, Path] = {}
    undesignated: list[Path] = []

    for path in vcf_paths:
        chromosomes = infer_vcf_chromosomes(path)
        if not chromosomes:
            undesignated.append(path)
            continue
        for chrom in chromosomes:
            index[chrom] = path

    fallback: Path | None = None
    if len(undesignated) == 1:
        fallback = undesignated[0]
    elif len(undesignated) > 1:
        LOG.warning(
            "Multiple VCFs without chromosome hints in filename; using first as fallback: %s",
            undesignated[0],
        )
        fallback = undesignated[0]

    return index, fallback


def read_manifest_csv(path: Path) -> pd.DataFrame:
    """Load the SNP manifest CSV produced by :func:`build_manifest`."""
    df = pd.read_csv(path)
    required = {"SNP_rsID", "Chromosome", "Position_GRCh38"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Manifest CSV is missing required columns: {sorted(missing)}")
    return df


class VcfHandleCache:
    """Lazy cache of open :class:`cyvcf2.VCF` readers."""

    def __init__(self) -> None:
        self._open: dict[str, cyvcf2.VCF] = {}

    def get(self, path: Path) -> cyvcf2.VCF:
        key = str(path.resolve())
        if key not in self._open:
            self._open[key] = cyvcf2.VCF(key)
        return self._open[key]

    def close(self) -> None:
        for handle in self._open.values():
            handle.close()
        self._open.clear()


def fetch_variant_at_locus(vcf: cyvcf2.VCF, chromosome: str, position: int) -> cyvcf2.Variant | None:
    """Return the variant record at ``chromosome:position`` using tabix, if present."""
    for alias in chromosome_query_aliases(chromosome):
        region = f"{alias}:{position}-{position}"
        try:
            for record in vcf(region):
                if int(record.POS) == position:
                    return record
        except (OSError, RuntimeError, ValueError) as exc:
            LOG.debug("Region query failed for %s: %s", region, exc)
    return None


def compute_cohort_allele_frequency(record: cyvcf2.Variant) -> tuple[float | None, int]:
    """Compute alternate-allele frequency and number of called samples."""
    n_called = 0
    n_alt_alleles = 0
    for genotype in record.genotypes:
        allele1, allele2 = genotype[0], genotype[1]
        if allele1 < 0 or allele2 < 0:
            continue
        n_called += 1
        if allele1 > 0:
            n_alt_alleles += 1
        if allele2 > 0:
            n_alt_alleles += 1

    n_alleles = n_called * 2
    if n_alleles == 0:
        return None, 0
    return n_alt_alleles / n_alleles, n_called


def resolve_vcf_for_chromosome(
    chromosome: str,
    *,
    chromosome_index: dict[str, Path],
    fallback_vcf: Path | None,
) -> Path | None:
    """Pick the VCF path that should contain variants on ``chromosome``."""
    normalised = normalize_chromosome_label(chromosome)
    if normalised in chromosome_index:
        return chromosome_index[normalised]
    return fallback_vcf


def extract_1kg_frequencies(manifest: pd.DataFrame, vcf_paths: Sequence[Path]) -> pd.DataFrame:
    """Extract 1KG allele frequencies for manifest SNPs via indexed region queries."""
    chromosome_index, fallback_vcf = build_chromosome_vcf_index(vcf_paths)
    cache = VcfHandleCache()

    rows: list[dict[str, Any]] = []
    try:
        for _, entry in manifest.iterrows():
            rs_id = str(entry["SNP_rsID"]).strip() if pd.notna(entry["SNP_rsID"]) else ""
            chrom_raw = entry["Chromosome"]
            pos_raw = entry["Position_GRCh38"]

            if not rs_id:
                LOG.warning("Skipping manifest row with empty rsID.")
                continue

            if pd.isna(chrom_raw) or pd.isna(pos_raw) or str(chrom_raw).strip() == "":
                LOG.warning("Manifest row %s lacks GRCh38 coordinates; skipping VCF lookup.", rs_id)
                rows.append(
                    {
                        "rsID": rs_id,
                        "chrom": "",
                        "pos": None,
                        "ref": "",
                        "alt": "",
                        "AF": None,
                        "N_called": 0,
                    }
                )
                continue

            chrom = normalize_chromosome_label(str(chrom_raw))
            position = int(float(pos_raw))
            vcf_path = resolve_vcf_for_chromosome(
                chrom,
                chromosome_index=chromosome_index,
                fallback_vcf=fallback_vcf,
            )
            if vcf_path is None:
                LOG.warning("No 1KG VCF found for chromosome %s (rsID %s).", chrom, rs_id)
                rows.append(
                    {
                        "rsID": rs_id,
                        "chrom": chrom,
                        "pos": int(position),
                        "ref": "",
                        "alt": "",
                        "AF": None,
                        "N_called": 0,
                    }
                )
                continue

            vcf = cache.get(vcf_path)
            record = fetch_variant_at_locus(vcf, chrom, position)
            if record is None:
                LOG.warning(
                    "SNP %s not found at %s:%s in %s.",
                    rs_id,
                    chrom,
                    position,
                    vcf_path.name,
                )
                rows.append(
                    {
                        "rsID": rs_id,
                        "chrom": chrom,
                        "pos": int(position),
                        "ref": "",
                        "alt": "",
                        "AF": None,
                        "N_called": 0,
                    }
                )
                continue

            ref = str(record.REF)
            alt = str(record.ALT[0]) if record.ALT else ""
            af, n_called = compute_cohort_allele_frequency(record)
            rows.append(
                {
                    "rsID": rs_id,
                    "chrom": chrom,
                    "pos": int(position),
                    "ref": ref,
                    "alt": alt,
                    "AF": af,
                    "N_called": n_called,
                }
            )
    finally:
        cache.close()

    return pd.DataFrame(rows)


def _add_build_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("overlapping_genes_with_snps.xlsx"),
        help="Path to overlapping_genes_with_snps.xlsx",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_MANIFEST_CSV,
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


def _add_extract_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST_CSV,
        help="Manifest CSV from the build subcommand",
    )
    parser.add_argument(
        "--vcf-glob",
        required=True,
        help="Path or glob to 1000 Genomes GRCh38 VCFs (bgzipped + tabix-indexed)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_1KG_FREQ_CSV,
        help="Output per-SNP allele-frequency table",
    )


def _add_log_level_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging verbosity",
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    argv = list(sys.argv[1:] if argv is None else argv)

    root = argparse.ArgumentParser(
        description=(
            "Build a GRCh38 SNP manifest from Ensembl (build) or extract 1000 "
            "Genomes allele frequencies for manifest SNPs (extract)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = root.add_subparsers(dest="command")

    build_parser = subparsers.add_parser(
        "build",
        help="Resolve rsIDs via Ensembl and write a GRCh38 manifest CSV",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_build_arguments(build_parser)
    _add_log_level_argument(build_parser)

    extract_parser = subparsers.add_parser(
        "extract",
        help="Extract 1KG allele frequencies for manifest SNPs from indexed VCFs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_extract_arguments(extract_parser)
    _add_log_level_argument(extract_parser)

    if argv and argv[0] in {"build", "extract"}:
        return root.parse_args(argv)

    if not argv or argv[0] in {"-h", "--help"}:
        root.print_help()
        raise SystemExit(0)

    args = build_parser.parse_args(argv)
    args.command = "build"
    return args


def run_build(args: argparse.Namespace) -> int:
    """Execute the manifest build workflow."""
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


def run_extract(args: argparse.Namespace) -> int:
    """Execute the 1KG allele-frequency extraction workflow."""
    manifest_path: Path = args.manifest
    output_path: Path = args.output

    if not manifest_path.is_file():
        LOG.error("Manifest CSV does not exist: %s", manifest_path.resolve())
        return 1

    try:
        vcf_paths = expand_vcf_paths(str(args.vcf_glob))
    except FileNotFoundError as exc:
        LOG.error("%s", exc)
        return 1

    if not vcf_paths:
        LOG.error("No VCF files matched: %s", args.vcf_glob)
        return 1

    LOG.info("Reading manifest from %s", manifest_path.resolve())
    manifest = read_manifest_csv(manifest_path)
    LOG.info(
        "Loaded %s manifest rows; querying %s VCF file(s)",
        len(manifest),
        len(vcf_paths),
    )

    frequencies = extract_1kg_frequencies(manifest, vcf_paths)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frequencies.to_csv(output_path, index=False)

    n_found = int(frequencies["AF"].notna().sum()) if not frequencies.empty else 0
    LOG.info(
        "Wrote %s SNP rows (%s with AF) to %s",
        len(frequencies),
        n_found,
        output_path.resolve(),
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    """CLI entry point: build manifest CSV or extract 1KG allele frequencies."""
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, str(args.log_level)),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    if args.command == "extract":
        return run_extract(args)
    return run_build(args)


if __name__ == "__main__":
    sys.exit(main())
