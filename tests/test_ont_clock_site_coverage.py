"""ONT versus bisulfite clock-CpG comparison and Zhou coordinate locks."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import polars as pl
import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SCRIPT = _REPO_ROOT / "scripts" / "clock" / "ont_clock_site_coverage.py"
_FIVE = _REPO_ROOT / "tests" / "fixtures" / "hm450_hg38_five_chr20.tsv"
_ZHOU = _REPO_ROOT / "data" / "dryrun" / "reference" / "HM450.hg38.manifest.tsv.gz"


def _load_script() -> ModuleType:
    name = "ont_clock_site_coverage"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


cov = _load_script()

_GRCH38_HEADER = "\n".join(
    [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr20\tLN:64444167\tUR:GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
        "@RG\tID:run_dna_r10.4.1_e8.2_400bps_sup@v5.0.0\tSM:hg002\tPL:ONT",
        "@PG\tID:basecaller\tPN:dorado",
        "",
    ]
)


def test_header_accepts_grch38_chr20() -> None:
    cov.require_grch38_chr20(_GRCH38_HEADER)


def test_header_stops_when_contig_is_not_chr20() -> None:
    header = _GRCH38_HEADER.replace("SN:chr20", "SN:20")
    with pytest.raises(ValueError, match="not named chr20"):
        cov.require_grch38_chr20(header)


def test_header_stops_when_reference_is_not_grch38() -> None:
    header = _GRCH38_HEADER.replace("GRCh38_no_alt", "hs37d5").replace("GCA_000001405.15_", "")
    with pytest.raises(ValueError, match="not GRCh38"):
        cov.require_grch38_chr20(header)


def test_bisulfite_bare_contig_20_is_refused() -> None:
    with pytest.raises(ValueError, match="not chr20"):
        cov.require_chr20_contig_name({"1", "20"})


def test_bismark_merges_plus_at_p_with_minus_at_p_plus_one() -> None:
    coverage = pl.DataFrame(
        {
            "chrom": ["chr20", "chr20", "chr20"],
            "start": [100, 101, 200],
            "n_meth": [2, 4, 1],
            "n_unmeth": [1, 1, 1],
        }
    )
    merged = cov.merge_bismark_cpg_strands(coverage).sort("start")
    assert merged.height == 2
    paired = merged.row(0, named=True)
    assert paired["start"] == 99
    assert paired["n_meth"] == 6
    assert paired["n_unmeth"] == 2
    assert paired["fraction"] == pytest.approx(0.75)
    singleton = merged.row(1, named=True)
    assert singleton["start"] == 199


def test_ont_fraction_adds_5mc_and_5hmc() -> None:
    pileup = pl.DataFrame(
        {
            "chrom": ["chr20", "chr20"],
            "start": [99, 200],
            "n_valid": [10, 0],
            "n_5mc": [3, 1],
            "n_5hmc": [2, 1],
        }
    )
    fractions = cov.ont_modification_fraction(pileup)
    assert fractions.height == 1
    assert fractions["fraction"][0] == pytest.approx(0.5)


def test_public_bedmethyl_sum_does_not_replace_own_pileup() -> None:
    own = pl.DataFrame(
        {
            "chrom": ["chr20"],
            "start": [99],
            "n_valid": [10],
            "n_5mc": [4],
            "n_5hmc": [1],
        }
    )
    before = own.clone()
    hap_1 = pl.DataFrame(
        {"chrom": ["chr20"], "start": [99], "n_valid": [6], "n_5mc": [2], "n_5hmc": [1]}
    )
    hap_2 = pl.DataFrame(
        {"chrom": ["chr20"], "start": [99], "n_valid": [3], "n_5mc": [2], "n_5hmc": [0]}
    )
    ungrouped = pl.DataFrame(
        {"chrom": ["chr20"], "start": [99], "n_valid": [1], "n_5mc": [0], "n_5hmc": [0]}
    )
    public = cov.sum_haplotype_bedmethyl([hap_1, hap_2, ungrouped])
    report = cov.pileup_agreement(own, public)
    assert own.equals(before)
    assert report.n_sites == 1
    assert report.n_exact == 1
    assert report.max_abs_count_delta == 0


def test_clock_cpg_comparison_reports_r_and_mean_absolute_difference() -> None:
    ont = pl.DataFrame(
        {
            "chrom": ["chr20", "chr20", "chr20"],
            "start": [99, 149, 199],
            "fraction": [0.50, 0.00, 1.00],
        }
    )
    bisulfite = pl.DataFrame(
        {
            "chrom": ["chr20", "chr20", "chr20"],
            "start": [99, 149, 199],
            "fraction": [0.40, 0.10, 0.90],
        }
    )
    clock = pl.DataFrame(
        {
            "chrom": ["chr20", "chr20", "chr20", "chr20"],
            "start": [99, 149, 199, 249],
            "probe_id": ["cgA", "cgB", "cgC", "cgMissing"],
        }
    )
    result = cov.compare_clock_cpgs(ont, bisulfite, clock)
    assert result.n_compared == 3
    assert result.pearson_r == pytest.approx(0.9897433186)
    assert result.mean_absolute_difference == pytest.approx(0.1)


def test_downsample_spec_records_seed_zero() -> None:
    assert cov.ILLUMINA_DOWNSAMPLE_SPEC == "0.1"
    assert cov.downsample_seed(cov.ILLUMINA_DOWNSAMPLE_SPEC) == 0
    assert cov.ILLUMINA_DOWNSAMPLE_SEED == 0
    assert "extra_analysis/all.bam" in cov.FORBIDDEN_BAM


def test_five_known_probes_match_modkit_starts() -> None:
    manifest = cov.load_zhou_manifest(_FIVE)
    cov.assert_known_zhou_probes(manifest)
    for probe_id, _chrom, cpg_beg in cov.KNOWN_CHR20_PROBES:
        assert cov.zhou_modkit_start(cpg_beg) == cpg_beg
        assert probe_id.startswith("cg")


@pytest.mark.skipif(not _ZHOU.is_file(), reason="Zhou manifest download is not in data/dryrun")
def test_five_known_probes_match_the_downloaded_zhou_manifest() -> None:
    manifest = cov.load_zhou_manifest(_ZHOU)
    cov.assert_known_zhou_probes(manifest)
