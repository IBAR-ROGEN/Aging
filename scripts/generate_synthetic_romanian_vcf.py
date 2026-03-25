#!/usr/bin/env python3
"""Stream a synthetic VCF v4.2 for a mock Romanian (EUR-structured) cohort.

Writes variant records one line at a time (no full-file buffer). Genotypes are
drawn under Hardy–Weinberg equilibrium using per-site alternate-allele
frequencies in [0.01, 0.5], mimicking common European allele-frequency ranges.

Example:
    uv run scripts/generate_synthetic_romanian_vcf.py --samples 50 --variants 200 \\
        --output /tmp/mock_ro.vcf
    bcftools view -H /tmp/mock_ro.vcf | head
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import TextIO

import numpy as np

# GRCh38 primary assembly lengths (bp) for chr1–chr22 — used in ##contig headers.
GRCH38_CHROM_LENGTHS: dict[str, int] = {
    "chr1": 248_956_422,
    "chr2": 242_193_529,
    "chr3": 198_295_559,
    "chr4": 190_214_555,
    "chr5": 181_538_259,
    "chr6": 170_805_979,
    "chr7": 159_345_973,
    "chr8": 145_138_636,
    "chr9": 138_394_717,
    "chr10": 133_797_422,
    "chr11": 135_086_622,
    "chr12": 133_275_309,
    "chr13": 114_364_328,
    "chr14": 107_043_718,
    "chr15": 101_991_189,
    "chr16": 90_338_345,
    "chr17": 83_257_508,
    "chr18": 80_373_285,
    "chr19": 58_617_616,
    "chr20": 64_444_167,
    "chr21": 46_709_983,
    "chr22": 50_818_468,
}

CHROMOSOMES: tuple[str, ...] = tuple(GRCH38_CHROM_LENGTHS.keys())


def configure_logging(level: int = logging.INFO) -> None:
    """Configure the root logger for stderr output in a standard format.

    Args:
        level: Logging level (e.g. ``logging.INFO``).
    """
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
        stream=sys.stderr,
    )


def hardy_weinberg_genotype_probs(alt_allele_freq: float) -> tuple[float, float, float]:
    """Return Hardy–Weinberg probabilities for a biallelic diploid site.

    Let ``p`` be the reference allele frequency and ``q = 1 - p`` the alternate
    frequency. Returns ``(P(REF/REF), P(REF/ALT), P(ALT/ALT))``.

    Args:
        alt_allele_freq: Frequency of the alternate allele ``q``, in ``(0, 0.5]``.

    Returns:
        Tuple of probabilities summing to 1.0 for hom-ref, het, and hom-alt.
    """
    q = float(alt_allele_freq)
    p = 1.0 - q
    hom_ref = p * p
    het = 2.0 * p * q
    hom_alt = q * q
    return hom_ref, het, hom_alt


def draw_genotypes_hardy_weinberg(
    rng: np.random.Generator,
    alt_allele_freq: float,
    n_samples: int,
) -> np.ndarray:
    """Draw diploid genotype codes under Hardy–Weinberg equilibrium.

    Args:
        rng: NumPy random generator instance.
        alt_allele_freq: Alternate allele frequency ``q``.
        n_samples: Number of individuals.

    Returns:
        One-dimensional array of dtype ``int64`` with values 0, 1, or 2.
    """
    hom_ref, het, hom_alt = hardy_weinberg_genotype_probs(alt_allele_freq)
    probs = np.array([hom_ref, het, hom_alt], dtype=np.float64)
    return rng.choice(np.array([0, 1, 2], dtype=np.int64), size=n_samples, p=probs)


def random_snp_alleles(rng: np.random.Generator) -> tuple[str, str]:
    """Sample a biallelic SNP reference and alternate base (both single uppercase).

    Args:
        rng: NumPy random generator instance.

    Returns:
        ``(ref, alt)`` with ``ref != alt``.
    """
    bases = ("A", "C", "G", "T")
    i, j = rng.choice(4, size=2, replace=False)
    return bases[int(i)], bases[int(j)]


def variant_chrom_and_pos(variant_index: int, n_variants: int) -> tuple[str, int]:
    """Assign chromosome and position in genome-sorted order (bcftools/tabix friendly).

    Variants are laid out in contiguous blocks per chromosome (chr1, then chr2,
    …) with strictly increasing positions within each chromosome so the emitted
    file is sorted as required by ``bcftools index`` without an extra sort pass.

    Args:
        variant_index: Zero-based variant index.
        n_variants: Total variant count (used to spread positions safely).

    Returns:
        ``(chrom, pos)`` with ``1 <= pos <= chrom_length``.
    """
    n_chrom = len(CHROMOSOMES)
    per_chrom = max(1, (n_variants + n_chrom - 1) // n_chrom)
    chrom_idx = min(variant_index // per_chrom, n_chrom - 1)
    local_index = variant_index - chrom_idx * per_chrom
    chrom = CHROMOSOMES[chrom_idx]
    chrom_len = GRCH38_CHROM_LENGTHS[chrom]
    stride = max(50, chrom_len // max(per_chrom, 1))
    pos = 10_000 + local_index * stride
    pos = int(min(max(pos, 1), chrom_len))
    return chrom, pos


def simulate_ad_dp_gq(
    rng: np.random.Generator,
    genotype: int,
    mean_depth: float,
) -> tuple[int, int, int, int]:
    """Simulate allelic depths, total depth, and genotype quality for one sample.

    Depth is Poisson-distributed; allelic counts follow binomial models
    consistent with the called genotype. ``GQ`` is a bounded integer suitable
    for FORMAT ``GQ``.

    Args:
        rng: NumPy random generator instance.
        genotype: ``0`` (hom ref), ``1`` (het), or ``2`` (hom alt).
        mean_depth: Mean total read depth.

    Returns:
        ``(ref_ad, alt_ad, dp, gq)``.
    """
    dp = int(max(1, rng.poisson(mean_depth)))
    error_rate = 0.005

    if genotype == 0:
        alt_ad = int(rng.binomial(dp, error_rate))
        ref_ad = dp - alt_ad
        gq = min(99, max(10, int(30 + 70 * (1.0 - error_rate))))
    elif genotype == 2:
        ref_ad = int(rng.binomial(dp, error_rate))
        alt_ad = dp - ref_ad
        gq = min(99, max(10, int(30 + 70 * (1.0 - error_rate))))
    else:
        ref_ad = int(rng.binomial(dp, 0.5))
        alt_ad = dp - ref_ad
        if ref_ad == 0 or alt_ad == 0:
            gq = min(99, max(15, dp))
        else:
            p_err = min(0.25, 2.0 * min(ref_ad, alt_ad) / max(dp, 1))
            gq = min(99, max(15, int(-10.0 * np.log10(max(p_err, 1e-9)))))

    return ref_ad, alt_ad, dp, gq


def format_gt_field(genotype: int) -> str:
    """Map genotype code to VCF ``GT`` string (unphased, diploid).

    Args:
        genotype: ``0``, ``1``, or ``2``.

    Returns:
        ``\"0/0\"``, ``\"0/1\"``, or ``\"1/1\"``.
    """
    if genotype == 0:
        return "0/0"
    if genotype == 1:
        return "0/1"
    return "1/1"


def format_sample_column(
    rng: np.random.Generator,
    genotype: int,
    mean_depth: float,
) -> str:
    """Build one sample's FORMAT values ``GT:AD:DP:GQ``.

    Args:
        rng: NumPy random generator instance.
        genotype: Genotype code 0, 1, or 2.
        mean_depth: Mean read depth for simulation.

    Returns:
        Single FORMAT field string for one sample.
    """
    ref_ad, alt_ad, dp, gq = simulate_ad_dp_gq(rng, genotype, mean_depth)
    gt = format_gt_field(genotype)
    return f"{gt}:{ref_ad},{alt_ad}:{dp}:{gq}"


def allele_counts_from_genotypes(genotypes: np.ndarray) -> tuple[int, int]:
    """Compute alternate allele count and total allele number from genotype codes.

    Encoding: ``0`` = two ref alleles, ``1`` = one ref and one alt, ``2`` = two alt.

    Args:
        genotypes: One-dimensional array with values 0, 1, or 2.

    Returns:
        ``(AC, AN)`` with ``AN = 2 * n_samples`` and ``AC`` the total alt count.
    """
    ac = int(np.sum(genotypes, dtype=np.int64))
    an = int(2 * genotypes.size)
    return ac, an


def write_vcf_headers(
    out: TextIO,
    sample_ids: Sequence[str],
    cohort_label: str,
) -> None:
    """Write VCF v4.2 meta and column header lines.

    Args:
        out: Output text stream.
        sample_ids: Sample column names in order.
        cohort_label: Short label for ``##cohort`` meta line.
    """
    out.write("##fileformat=VCFv4.2\n")
    out.write(f"##source=generate_synthetic_romanian_vcf.py\n")
    out.write(f"##synthetic_cohort={cohort_label}\n")
    out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    out.write(
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">\n'
    )
    out.write(
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">\n'
    )
    out.write(
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes">\n'
    )
    out.write(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n'
    )
    out.write(
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    )
    out.write(
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n'
    )
    out.write(
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">\n'
    )
    out.write(
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n'
    )
    for chrom, length in GRCH38_CHROM_LENGTHS.items():
        out.write(f"##contig=<ID={chrom},length={length}>\n")

    tabs = "\t".join(sample_ids)
    out.write(
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{tabs}\n"
    )


def romanian_cohort_sample_ids(n_samples: int, prefix: str = "RO_EUR") -> list[str]:
    """Build mock Romanian cohort sample names (bcftools-safe identifiers).

    Args:
        n_samples: Number of samples.
        prefix: Alphanumeric prefix for each ID.

    Returns:
        List of length ``n_samples`` of distinct sample column names.
    """
    width = max(6, len(str(n_samples)))
    return [f"{prefix}_{i:0{width}d}" for i in range(1, n_samples + 1)]


def draw_site_alt_frequency(rng: np.random.Generator) -> float:
    """Draw a synthetic European-like alternate allele frequency for one site.

    Samples uniformly from ``[0.01, 0.5]``, the range requested for EUR-like
    common-variant simulation.

    Args:
        rng: NumPy random generator instance.

    Returns:
        Alternate allele frequency ``q`` for Hardy–Weinberg sampling.
    """
    return float(rng.uniform(0.01, 0.5))


def iter_variant_lines(
    rng: np.random.Generator,
    n_samples: int,
    n_variants: int,
    mean_depth: float,
) -> Iterator[str]:
    """Yield VCF data lines (no header), one variant per iteration.

    Memory use is ``O(n_samples)`` per variant for genotype and FORMAT columns,
    not ``O(n_variants * n_samples)`` for the full file.

    Args:
        rng: NumPy random generator instance.
        n_samples: Number of diploid samples.
        n_variants: Number of variant rows to emit.
        mean_depth: Mean read depth for ``AD``/``DP`` simulation.

    Yields:
        Single tab-separated VCF variant line including trailing newline.
    """
    log = logging.getLogger(__name__)
    fmt = "GT:AD:DP:GQ"
    for variant_index in range(n_variants):
        chrom, pos = variant_chrom_and_pos(variant_index, n_variants)
        ref, alt = random_snp_alleles(rng)
        alt_freq = draw_site_alt_frequency(rng)
        genotypes = draw_genotypes_hardy_weinberg(rng, alt_freq, n_samples)
        ac, an = allele_counts_from_genotypes(genotypes)
        af_info = (ac / an) if an else 0.0
        info = f"AC={ac};AN={an};AF={af_info:.6f};END={pos}"
        parts: list[str] = [
            chrom,
            str(pos),
            f"RO_MOCK_{variant_index + 1:08d}",
            ref,
            alt,
            "60",
            "PASS",
            info,
            fmt,
        ]
        for sample_idx in range(n_samples):
            parts.append(
                format_sample_column(
                    rng, int(genotypes[sample_idx]), mean_depth
                )
            )
        line = "\t".join(parts) + "\n"
        if variant_index > 0 and variant_index % 10_000 == 0:
            log.info("Emitted %d variant lines", variant_index)
        yield line


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for synthetic VCF generation.

    Args:
        argv: Argument list; defaults to ``sys.argv[1:]`` when ``None``.

    Returns:
        Parsed namespace with ``samples``, ``variants``, ``output``, etc.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Generate a synthetic VCF v4.2 for a mock Romanian (EUR-structured) "
            "population cohort under Hardy–Weinberg genotypes."
        )
    )
    parser.add_argument(
        "--samples",
        type=int,
        required=True,
        metavar="N",
        help="Number of diploid samples (columns after FORMAT).",
    )
    parser.add_argument(
        "--variants",
        type=int,
        required=True,
        metavar="M",
        help="Number of variant rows to write.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output path for uncompressed VCF (tab-delimited text).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility.",
    )
    parser.add_argument(
        "--mean-depth",
        type=float,
        default=32.0,
        help="Mean simulated read depth per sample per site (Poisson mean).",
    )
    parser.add_argument(
        "--cohort-label",
        type=str,
        default="mock_RO_EUR_cohort",
        help="Label written to ##synthetic_cohort header line.",
    )
    parser.add_argument(
        "--sample-prefix",
        type=str,
        default="RO_EUR",
        help="Prefix for synthetic sample IDs.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    """CLI entry: write headers then stream variant lines to the output file.

    Args:
        argv: Optional argument list for testing; defaults to process arguments.
    """
    args = parse_args(argv)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    configure_logging(log_level)
    log = logging.getLogger(__name__)

    if args.samples < 1:
        raise SystemExit("--samples must be >= 1")
    if args.variants < 1:
        raise SystemExit("--variants must be >= 1")
    if args.mean_depth <= 0:
        raise SystemExit("--mean-depth must be positive")

    rng = np.random.default_rng(args.seed)
    sample_ids = romanian_cohort_sample_ids(args.samples, args.sample_prefix)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    log.info(
        "Writing %d variants for %d samples to %s",
        args.variants,
        args.samples,
        args.output,
    )
    with args.output.open("w", encoding="ascii", newline="\n") as out:
        write_vcf_headers(out, sample_ids, args.cohort_label)
        for line in iter_variant_lines(
            rng, args.samples, args.variants, args.mean_depth
        ):
            out.write(line)
    log.info("Finished writing VCF.")


if __name__ == "__main__":
    main()
