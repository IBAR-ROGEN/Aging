"""Synthetic UKB-RAP folder generator (phenotypes + LA-SNP VCF)."""

from __future__ import annotations

import logging
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl
import typer

from rogen_aging.vcf import (
    GRCH38_CHROM_LENGTHS,
    allele_counts_from_genotypes,
    configure_logging,
    draw_genotypes_hardy_weinberg,
    draw_site_alt_frequency,
    format_sample_column,
    random_snp_alleles,
    write_vcf_headers,
)

ACTIVITY_HEADER = (
    "Activity 2.1.8.1 synthetic UKB-RAP mock cohort - safe for GitHub (no real EIDs)"
)

# v2 phenotype dictionary fields (excluding join key ``eid``).
PHENOTYPE_V2_FIELDS: tuple[str, ...] = (
    "age",
    "sex",
    "parental_longevity",
    "ad_diagnosis_code",
    "pd_diagnosis_code",
    "frailty_weight_loss",
    "frailty_exhaustion",
    "frailty_weakness",
    "frailty_slowness",
    "frailty_low_activity",
)

DEFAULT_MANIFEST = Path("analysis/ukb_snp_manifest_v0.1.csv")
DEFAULT_OUTPUT_DIR = Path("test_data/mock_ukb_rap")
PHENOTYPE_REL_PATH = Path("phenotypes/ukb_phenotypes.csv")
GENOTYPE_REL_PATH = Path("genotypes/ukb_la_snps.vcf")

AD_ICD10_CODES: tuple[str, ...] = ("G30.0", "G30.1", "G30.8", "G30.9", "F00.0", "F00.1")
PD_ICD10_CODES: tuple[str, ...] = ("G20", "G21.0", "G21.1", "G21.2", "G22")

app = typer.Typer(
    help="Generate synthetic UKB-RAP-style phenotype table and LA-SNP VCF (Activity 2.1.8.1)."
)


@dataclass(frozen=True)
class ManifestSnp:
    """One LA-SNP row from the UKB SNP manifest."""

    rs_id: str
    chromosome: str
    position: int


def normalize_chromosome(chrom: str) -> str:
    """Map manifest chromosome labels to VCF ``#CHROM`` names (``chr1`` … ``chr22``)."""
    raw = chrom.strip()
    if raw.lower().startswith("chr"):
        return raw if raw.startswith("chr") else f"chr{raw[3:]}"
    return f"chr{raw}"


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Sort key for GRCh38 autosomes chr1–chr22."""
    normalized = normalize_chromosome(chrom)
    if normalized.startswith("chr"):
        body = normalized[3:]
        if body.isdigit():
            return int(body), normalized
    return 99, normalized


def load_snp_manifest(path: Path) -> list[ManifestSnp]:
    """Load and validate the UKB LA-SNP manifest CSV."""
    if not path.is_file():
        raise FileNotFoundError(f"SNP manifest not found: {path.resolve()}")

    frame = pl.read_csv(path)
    required = {"SNP_rsID", "Chromosome", "Position_GRCh38"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"SNP manifest missing columns: {sorted(missing)}")

    rows: list[ManifestSnp] = []
    for record in frame.iter_rows(named=True):
        rs_id = str(record["SNP_rsID"]).strip()
        chrom = str(record["Chromosome"]).strip()
        pos_raw = record["Position_GRCh38"]
        if not rs_id or not chrom or pos_raw is None:
            continue
        position = int(pos_raw)
        rows.append(ManifestSnp(rs_id=rs_id, chromosome=chrom, position=position))

    if not rows:
        raise ValueError(f"SNP manifest contains no usable rows: {path}")

    rows.sort(key=lambda snp: (chrom_sort_key(snp.chromosome), snp.position, snp.rs_id))
    return rows


def synthetic_eids(n_samples: int, *, prefix: str = "SYN_EID") -> list[str]:
    """Build synthetic participant IDs (not real UKB 7-digit EIDs)."""
    width = max(8, len(str(n_samples)))
    return [f"{prefix}_{i:0{width}d}" for i in range(1, n_samples + 1)]


def generate_phenotype_table(
    rng: np.random.Generator,
    eids: Sequence[str],
    *,
    min_age: int = 40,
    max_age: int = 80,
) -> pl.DataFrame:
    """Simulate the v2 phenotype dictionary for ``eids``."""
    n_samples = len(eids)
    age = rng.integers(min_age, max_age + 1, size=n_samples)
    sex = rng.integers(0, 2, size=n_samples)
    parental_longevity = rng.binomial(1, 0.15, size=n_samples)

    ad_flag = rng.binomial(1, 0.02, size=n_samples)
    pd_flag = rng.binomial(1, 0.008, size=n_samples)
    ad_codes = [
        rng.choice(AD_ICD10_CODES) if flag else "" for flag in ad_flag
    ]
    pd_codes = [
        rng.choice(PD_ICD10_CODES) if flag else "" for flag in pd_flag
    ]

    frailty_probs = rng.uniform(0.05, 0.25, size=5)
    frailty_cols = {
        field: rng.binomial(1, prob, size=n_samples)
        for field, prob in zip(
            (
                "frailty_weight_loss",
                "frailty_exhaustion",
                "frailty_weakness",
                "frailty_slowness",
                "frailty_low_activity",
            ),
            frailty_probs,
            strict=True,
        )
    }

    return pl.DataFrame(
        {
            "eid": list(eids),
            "age": age,
            "sex": sex,
            "parental_longevity": parental_longevity,
            "ad_diagnosis_code": ad_codes,
            "pd_diagnosis_code": pd_codes,
            **frailty_cols,
        }
    )


def write_phenotype_csv(path: Path, table: pl.DataFrame) -> None:
    """Write phenotype table with an Activity 2.1.8.1 safety header comment."""
    path.parent.mkdir(parents=True, exist_ok=True)
    csv_body = table.write_csv()
    path.write_text(f"# {ACTIVITY_HEADER}\n{csv_body}", encoding="utf-8")


def iter_manifest_variant_lines(
    rng: np.random.Generator,
    manifest: Sequence[ManifestSnp],
    sample_ids: Sequence[str],
    mean_depth: float,
) -> Iterator[str]:
    """Yield VCF variant lines for manifest SNPs (one line per SNP)."""
    n_samples = len(sample_ids)
    fmt = "GT:AD:DP:GQ"
    for snp in manifest:
        chrom = normalize_chromosome(snp.chromosome)
        if chrom not in GRCH38_CHROM_LENGTHS:
            raise ValueError(f"Unsupported chromosome in manifest: {snp.chromosome!r}")
        pos = snp.position
        ref, alt = random_snp_alleles(rng)
        alt_freq = draw_site_alt_frequency(rng)
        genotypes = draw_genotypes_hardy_weinberg(rng, alt_freq, n_samples)
        ac, an = allele_counts_from_genotypes(genotypes)
        af_info = (ac / an) if an else 0.0
        info = f"AC={ac};AN={an};AF={af_info:.6f};END={pos}"
        parts: list[str] = [
            chrom,
            str(pos),
            snp.rs_id,
            ref,
            alt,
            "60",
            "PASS",
            info,
            fmt,
        ]
        for sample_idx in range(n_samples):
            parts.append(
                format_sample_column(rng, int(genotypes[sample_idx]), mean_depth)
            )
        yield "\t".join(parts) + "\n"


def write_la_snp_vcf(
    path: Path,
    *,
    sample_ids: Sequence[str],
    manifest: Sequence[ManifestSnp],
    rng: np.random.Generator,
    mean_depth: float,
    cohort_label: str,
) -> None:
    """Stream LA-SNP VCF with headers and manifest-ordered variant rows."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as out:
        write_vcf_headers(
            out,
            sample_ids,
            cohort_label,
            source="rogen_aging.ukb.mock_rap",
            extra_meta_lines=(
                f"##activity={ACTIVITY_HEADER}",
                "##synthetic_data=true",
            ),
        )
        for line in iter_manifest_variant_lines(
            rng, manifest, sample_ids, mean_depth
        ):
            out.write(line)


def generate_ukb_rap_mock(
    *,
    n_samples: int,
    snp_manifest: Path,
    output_dir: Path,
    seed: int | None,
    mean_depth: float = 32.0,
) -> tuple[Path, Path]:
    """Create phenotype CSV and matching LA-SNP VCF under ``output_dir``."""
    if n_samples < 1:
        raise ValueError("--n-samples must be >= 1")
    if mean_depth <= 0:
        raise ValueError("mean_depth must be positive")

    rng = np.random.default_rng(seed)
    manifest = load_snp_manifest(snp_manifest)
    eids = synthetic_eids(n_samples)
    phenotype_path = output_dir / PHENOTYPE_REL_PATH
    genotype_path = output_dir / GENOTYPE_REL_PATH

    phenotype_table = generate_phenotype_table(rng, eids)
    write_phenotype_csv(phenotype_path, phenotype_table)
    write_la_snp_vcf(
        genotype_path,
        sample_ids=eids,
        manifest=manifest,
        rng=rng,
        mean_depth=mean_depth,
        cohort_label="mock_ukb_rap_la_snp_cohort",
    )
    return phenotype_path, genotype_path


@app.command()
def main(
    n_samples: int = typer.Option(
        1000,
        "--n-samples",
        "-n",
        help="Number of synthetic participants.",
    ),
    snp_manifest: Path = typer.Option(
        DEFAULT_MANIFEST,
        "--snp-manifest",
        exists=False,
        path_type=Path,
        help="UKB LA-SNP manifest CSV (rsID + GRCh38 coordinates).",
    ),
    output_dir: Path = typer.Option(
        DEFAULT_OUTPUT_DIR,
        "--output-dir",
        "-o",
        path_type=Path,
        help="Root directory for UKB-RAP-style layout.",
    ),
    seed: int | None = typer.Option(
        None,
        "--seed",
        "-s",
        help="Random seed for reproducibility.",
    ),
    mean_depth: float = typer.Option(
        32.0,
        "--mean-depth",
        help="Mean simulated read depth per sample per site (Poisson mean).",
    ),
    verbose: bool = typer.Option(
        False,
        "-v",
        "--verbose",
        help="Enable debug logging.",
    ),
) -> None:
    """Generate synthetic UKB-RAP phenotype table and LA-SNP VCF."""
    configure_logging(logging.DEBUG if verbose else logging.INFO)
    log = logging.getLogger(__name__)

    phenotype_path, genotype_path = generate_ukb_rap_mock(
        n_samples=n_samples,
        snp_manifest=snp_manifest,
        output_dir=output_dir,
        seed=seed,
        mean_depth=mean_depth,
    )
    manifest = load_snp_manifest(snp_manifest)
    log.info(
        "Wrote %d samples, %d LA-SNPs → %s and %s",
        n_samples,
        len(manifest),
        phenotype_path,
        genotype_path,
    )
    typer.echo(f"Phenotypes: {phenotype_path}")
    typer.echo(f"Genotypes:  {genotype_path}")


if __name__ == "__main__":
    app()
