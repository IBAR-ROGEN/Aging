#!/usr/bin/env python3
"""ONT versus bisulfite methylation at covered chr20 clock CpGs.

URL decisions. Do not use ``gm24385_mod_2021.09/extra_analysis/all.bam``.

* ONT reads: ``s3://ont-open-data/giab_2025.01/analysis/wf-human-variation/sup/HG002/PAW70337/output/SAMPLE.haplotagged.cram``
  (and ``.crai``). Header ``@SQ`` uses ``chr20`` on ``GCA_000001405.15_GRCh38_no_alt_analysis_set``.
  ``@RG`` records ``basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0`` and
  ``modbase_models=dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mC_5hmC@v2.0.1``,
  flow cell chemistry R10.4.1 / e8.2, instrument ``PC48LILITH``, run ``PAW70337``.
* Bisulfite, same cell line: ``s3://ont-open-data/gm24385_mod_2021.09/bisulphite/cpg/CpG.gz.bismark.cov.gz``.
  Contig names must be ``chr20`` before a join. ``.cov`` positions are 1-based and
  per strand: C(+) at ``p`` and C(-) at ``p+1`` become one CpG.
* Illumina: ``HG002.GRCh38.300x_chr20.bam`` in
  ``NHGRI_Illumina300X_AJtrio_novoalign_bams``. Downsample with
  ``samtools view -s 0.1`` (seed 0, fraction 0.1, about 30x) before writing FASTQs.
* Zhou HM450 hg38 manifest (``CpG_beg`` is 0-based and equals the modkit start):
  ``https://zhouserver.research.chop.edu/InfiniumAnnotation/current/HM450/HM450.hg38.manifest.tsv.gz``.

HG002, NA24385, and GM24385 are one cell line. The ONT run is dated 2024-06-19
and the bisulfite file was published in the 2021.09 bucket. The samples were
cultured and sequenced years apart, so part of any difference is biological.

Own modkit pileup stays the ONT measurement. Summing the three haplotype
bedMethyl files (``.1``, ``.2``, ``.ungrouped``) is a cross-check only.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl
import typer
from scipy.stats import pearsonr

ONT_CRAM_URL = (
    "https://ont-open-data.s3.amazonaws.com/giab_2025.01/analysis/"
    "wf-human-variation/sup/HG002/PAW70337/output/SAMPLE.haplotagged.cram"
)
BISULFITE_COV_URL = (
    "https://ont-open-data.s3.amazonaws.com/gm24385_mod_2021.09/bisulphite/cpg/"
    "CpG.gz.bismark.cov.gz"
)
ZHOU_HM450_HG38_URL = (
    "https://zhouserver.research.chop.edu/InfiniumAnnotation/current/HM450/"
    "HM450.hg38.manifest.tsv.gz"
)
ILLUMINA_CHR20_BAM_URL = (
    "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/"
    "HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/"
    "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x_chr20.bam"
)
# samtools view -s INT.FRAC: INT is the seed, FRAC is the keep fraction.
ILLUMINA_DOWNSAMPLE_SPEC = "0.1"
ILLUMINA_DOWNSAMPLE_SEED = 0
GRCH38_CHR20_LENGTH = 64_444_167
FORBIDDEN_BAM = "gm24385_mod_2021.09/extra_analysis/all.bam"

# Copied from HM450.hg38.manifest.tsv.gz. CpG_beg is the 0-based modkit start.
KNOWN_CHR20_PROBES: tuple[tuple[str, str, int], ...] = (
    ("cg14814714", "chr20", 86358),
    ("cg18865733", "chr20", 86440),
    ("cg08088390", "chr20", 87528),
    ("cg01958189", "chr20", 87754),
    ("cg12738399", "chr20", 92460),
)


@dataclass(frozen=True)
class CoverageComparison:
    """ONT fraction versus bisulfite fraction at covered clock CpGs."""

    n_compared: int
    pearson_r: float
    mean_absolute_difference: float


@dataclass(frozen=True)
class PileupAgreement:
    """Count-level agreement. Does not replace the local pileup."""

    n_sites: int
    n_exact: int
    max_abs_count_delta: int


def downsample_seed(spec: str) -> int:
    """Return the integer seed encoded by a samtools ``view -s`` value."""
    whole, dot, frac = spec.partition(".")
    if dot == "" or frac == "" or not frac.isdigit() or (whole != "" and not whole.isdigit()):
        raise ValueError(f"expected samtools -s seed.fraction, got {spec!r}")
    return int(whole) if whole else 0


def require_grch38_chr20(header_text: str) -> None:
    """Stop unless ``@SQ`` names ``chr20`` on GRCh38.

    Args:
        header_text: SAM header text (``@SQ``, ``@RG``, ``@PG``).

    Raises:
        ValueError: If chr20 is absent, not named ``chr20``, or the reference
            is not GRCh38.
    """
    sq_lines = [line for line in header_text.splitlines() if line.startswith("@SQ")]
    if not sq_lines:
        raise ValueError("CRAM header has no @SQ lines")
    chr20 = [line for line in sq_lines if "\tSN:chr20\t" in line]
    if len(chr20) != 1:
        raise ValueError("reference contig is not named chr20; refusing to continue")
    if f"LN:{GRCH38_CHR20_LENGTH}" not in chr20[0]:
        raise ValueError("chr20 length is not the GRCh38 length; refusing to continue")
    for line in sq_lines:
        if "UR:" not in line:
            continue
        if "GRCh38" not in line and "GCA_000001405.15" not in line:
            raise ValueError("reference is not GRCh38; refusing to continue")


def require_chr20_contig_name(chromosomes: set[str]) -> None:
    """Require the bisulfite file to call the contig ``chr20``.

    Args:
        chromosomes: Contig names observed in the coverage file.

    Raises:
        ValueError: If ``chr20`` is missing, including a bare ``20``.
    """
    if "chr20" in chromosomes:
        return
    if "20" in chromosomes:
        raise ValueError("bisulfite contig is 20, not chr20; refusing to join")
    raise ValueError("bisulfite file has no chr20 contig; refusing to join")


def merge_bismark_cpg_strands(coverage: pl.DataFrame) -> pl.DataFrame:
    """Merge per-strand Bismark ``.cov`` rows into one CpG.

    Positions are 1-based. C(+) at ``p`` and C(-) at ``p+1`` on the same contig
    are one CpG. The merged start is 0-based ``p - 1``, which matches modkit
    and Zhou ``CpG_beg``. Counts are summed. A row that is not one base after
    the previous unpaired row stays a singleton at ``start - 1``.

    Args:
        coverage: Columns ``chrom``, ``start``, ``n_meth``, ``n_unmeth``.

    Returns:
        One row per CpG with 0-based ``start`` and a fraction in ``[0, 1]``.
    """
    required = {"chrom", "start", "n_meth", "n_unmeth"}
    missing = required.difference(coverage.columns)
    if missing:
        raise ValueError(f"bismark coverage missing columns: {sorted(missing)}")
    ordered = coverage.sort(["chrom", "start"]).iter_rows(named=True)
    merged: list[dict[str, object]] = []
    pending: dict[str, object] | None = None
    for row in ordered:
        if pending is None:
            pending = row
            continue
        same_chrom = row["chrom"] == pending["chrom"]
        if same_chrom and int(row["start"]) == int(pending["start"]) + 1:
            merged.append(_merged_cpg(pending, row))
            pending = None
            continue
        merged.append(_merged_cpg(pending, None))
        pending = row
    if pending is not None:
        merged.append(_merged_cpg(pending, None))
    if not merged:
        return pl.DataFrame(
            schema={
                "chrom": pl.Utf8,
                "start": pl.Int64,
                "n_meth": pl.Int64,
                "n_unmeth": pl.Int64,
                "fraction": pl.Float64,
            }
        )
    return pl.DataFrame(merged)


def _merged_cpg(
    plus_or_single: dict[str, object],
    minus: dict[str, object] | None,
) -> dict[str, object]:
    n_meth = int(plus_or_single["n_meth"])
    n_unmeth = int(plus_or_single["n_unmeth"])
    if minus is not None:
        n_meth += int(minus["n_meth"])
        n_unmeth += int(minus["n_unmeth"])
    depth = n_meth + n_unmeth
    fraction = (n_meth / depth) if depth else float("nan")
    return {
        "chrom": str(plus_or_single["chrom"]),
        "start": int(plus_or_single["start"]) - 1,
        "n_meth": n_meth,
        "n_unmeth": n_unmeth,
        "fraction": fraction,
    }


def ont_modification_fraction(pileup: pl.DataFrame) -> pl.DataFrame:
    """Return ``(5mC + 5hmC) / valid`` at each 0-based start.

    Args:
        pileup: Columns ``chrom``, ``start``, ``n_valid``, ``n_5mc``, ``n_5hmc``.

    Returns:
        Rows with ``n_valid > 0`` and fraction ``(n_5mc + n_5hmc) / n_valid``.
    """
    required = {"chrom", "start", "n_valid", "n_5mc", "n_5hmc"}
    missing = required.difference(pileup.columns)
    if missing:
        raise ValueError(f"ONT pileup missing columns: {sorted(missing)}")
    covered = pileup.filter(pl.col("n_valid") > 0)
    return covered.with_columns(
        ((pl.col("n_5mc") + pl.col("n_5hmc")) / pl.col("n_valid")).alias("fraction")
    )


def sum_haplotype_bedmethyl(parts: list[pl.DataFrame]) -> pl.DataFrame:
    """Sum modkit counts across haplotype bedMethyl tables.

    The public files are ``SAMPLE.wf_mods.1``, ``.2``, and ``.ungrouped``.
    Fractions are not averaged. ``n_valid``, ``n_5mc``, and ``n_5hmc`` are summed.

    Args:
        parts: One frame per haplotype file, with those count columns.

    Returns:
        One row per ``chrom`` and 0-based ``start``.
    """
    if not parts:
        raise ValueError("no haplotype bedMethyl tables to sum")
    stacked = pl.concat(parts, how="vertical")
    return stacked.group_by(["chrom", "start"]).agg(
        pl.col("n_valid").sum(),
        pl.col("n_5mc").sum(),
        pl.col("n_5hmc").sum(),
    )


def pileup_agreement(own: pl.DataFrame, public_sum: pl.DataFrame) -> PileupAgreement:
    """Compare count columns. ``own`` is not replaced or rewritten.

    Args:
        own: Local modkit pileup counts.
        public_sum: Sum of the three ONT-provided haplotype bedMethyl files.

    Returns:
        How many shared sites match exactly, and the largest absolute count gap.
    """
    keys = ["chrom", "start"]
    counts = ["n_valid", "n_5mc", "n_5hmc"]
    joined = own.select([*keys, *counts]).join(
        public_sum.select([*keys, *counts]),
        on=keys,
        how="inner",
        suffix="_public",
    )
    if joined.height == 0:
        return PileupAgreement(n_sites=0, n_exact=0, max_abs_count_delta=0)
    delta = joined.select(
        [(pl.col(name) - pl.col(f"{name}_public")).abs().alias(name) for name in counts]
    )
    exact = delta.filter(
        (pl.col("n_valid") == 0) & (pl.col("n_5mc") == 0) & (pl.col("n_5hmc") == 0)
    )
    max_delta = max(int(delta[name].max()) for name in counts)
    return PileupAgreement(
        n_sites=joined.height,
        n_exact=exact.height,
        max_abs_count_delta=max_delta,
    )


def compare_clock_cpgs(
    ont: pl.DataFrame,
    bisulfite: pl.DataFrame,
    clock_cpgs: pl.DataFrame,
) -> CoverageComparison:
    """Pearson r and mean absolute difference at covered chr20 clock CpGs.

    Both fractions are on a 0–1 scale. Bisulfite percent methylation is
    ``n_meth / (n_meth + n_unmeth)`` after the strand merge, not the raw
    per-strand percentage. ONT is ``(5mC + 5hmC) / valid``.

    Args:
        ont: Output of :func:`ont_modification_fraction`.
        bisulfite: Output of :func:`merge_bismark_cpg_strands`, already limited
            to contigs that passed :func:`require_chr20_contig_name`.
        clock_cpgs: ``chrom``, ``start`` (0-based), and ``probe_id`` for clock CpGs.

    Returns:
        Correlation and mean absolute difference for sites covered by both.

    Raises:
        ValueError: If fewer than two clock CpGs are covered on both sides.
    """
    ont_frac = ont.select(["chrom", "start", pl.col("fraction").alias("ont")])
    bs_frac = bisulfite.select(["chrom", "start", pl.col("fraction").alias("bisulfite")])
    sites = clock_cpgs.select(["chrom", "start", "probe_id"])
    compared = sites.join(ont_frac, on=["chrom", "start"], how="inner").join(
        bs_frac, on=["chrom", "start"], how="inner"
    )
    compared = compared.filter(pl.col("ont").is_not_nan() & pl.col("bisulfite").is_not_nan())
    if compared.height < 2:
        raise ValueError(
            "fewer than two covered chr20 clock CpGs have both ONT and bisulfite calls"
        )
    ont_values = compared["ont"].to_list()
    bs_values = compared["bisulfite"].to_list()
    correlation = float(pearsonr(ont_values, bs_values).statistic)
    absolute = [abs(left - right) for left, right in zip(ont_values, bs_values, strict=True)]
    return CoverageComparison(
        n_compared=compared.height,
        pearson_r=correlation,
        mean_absolute_difference=sum(absolute) / len(absolute),
    )


def zhou_modkit_start(cpg_beg: int) -> int:
    """Return the modkit start for a Zhou ``CpG_beg``.

    ``CpG_beg`` is already 0-based. It is the start of the 2 bp CpG and the
    modkit coordinate of that cytosine. Do not add 1.
    """
    if cpg_beg < 0:
        raise ValueError(f"CpG_beg must be >= 0, got {cpg_beg}")
    return cpg_beg


def assert_known_zhou_probes(manifest: pl.DataFrame) -> None:
    """Check five chr20 probes whose ``CpG_beg`` equals the modkit start.

    Args:
        manifest: Zhou HM450 hg38 table with ``probeID``, ``CpG_chrm``, ``CpG_beg``.

    Raises:
        AssertionError: If a known probe is missing or its coordinate disagrees.
    """
    for probe_id, chrom, cpg_beg in KNOWN_CHR20_PROBES:
        rows = manifest.filter(pl.col("probeID") == probe_id)
        if rows.height != 1:
            raise AssertionError(f"{probe_id} should occur once, found {rows.height}")
        observed_chrom = str(rows["CpG_chrm"][0])
        observed_beg = int(rows["CpG_beg"][0])
        if observed_chrom != chrom or observed_beg != cpg_beg:
            raise AssertionError(
                f"{probe_id} is {observed_chrom}:{observed_beg}, expected {chrom}:{cpg_beg}"
            )
        if zhou_modkit_start(observed_beg) != cpg_beg:
            raise AssertionError(f"{probe_id} CpG_beg is not the modkit start")
        if int(rows["CpG_end"][0]) != cpg_beg + 2:
            raise AssertionError(f"{probe_id} CpG interval is not 2 bp")


def load_zhou_manifest(path: Path) -> pl.DataFrame:
    """Load the Zhou manifest TSV (plain or gzip)."""
    return pl.read_csv(
        path,
        separator="\t",
        columns=["CpG_chrm", "CpG_beg", "CpG_end", "probeID"],
        null_values=["NA"],
        schema_overrides={
            "CpG_chrm": pl.Utf8,
            "CpG_beg": pl.Int64,
            "CpG_end": pl.Int64,
            "probeID": pl.Utf8,
        },
    )


def main() -> None:
    """Print the recorded URL decisions and the Illumina downsample seed."""
    if FORBIDDEN_BAM in ONT_CRAM_URL or FORBIDDEN_BAM in BISULFITE_COV_URL:
        raise typer.Exit(code=2)
    typer.echo(f"ONT CRAM: {ONT_CRAM_URL}")
    typer.echo(f"Bisulfite cov: {BISULFITE_COV_URL}")
    typer.echo(f"Zhou manifest: {ZHOU_HM450_HG38_URL}")
    typer.echo(f"Illumina chr20 BAM: {ILLUMINA_CHR20_BAM_URL}")
    typer.echo(
        "Illumina downsample: samtools view -s "
        f"{ILLUMINA_DOWNSAMPLE_SPEC} (seed {downsample_seed(ILLUMINA_DOWNSAMPLE_SPEC)})"
    )
    typer.echo(f"Refusing: {FORBIDDEN_BAM}")


if __name__ == "__main__":
    typer.run(main)
