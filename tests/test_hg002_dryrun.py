"""Unit tests for HG002 short-read dry-run helpers."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SCRIPT = _REPO_ROOT / "scripts" / "dryrun" / "hg002_shortread_dryrun.py"
_SUMMARY_FIXTURE = _REPO_ROOT / "tests" / "fixtures" / "vcfeval_summary_small.txt"


def _load() -> ModuleType:
    name = "hg002_shortread_dryrun"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


dry = _load()


def test_illumina_bam_url_is_imported_not_copied() -> None:
    source = _SCRIPT.read_text(encoding="utf-8")
    assert "ILLUMINA_CHR20_BAM_URL" in source
    assert "ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio" not in source
    assert dry.ILLUMINA_BAM_URL.startswith("https://")
    assert "HG002.GRCh38.300x_chr20.bam" in dry.ILLUMINA_BAM_URL


def test_parse_region_full_and_window() -> None:
    full = dry.parse_region("chr20")
    assert full.chrom == "chr20"
    assert full.start is None and full.end is None
    assert full.as_samtools() == "chr20"

    window = dry.parse_region("chr20:10000000-20000000")
    assert window.chrom == "chr20"
    assert window.start == 10_000_000
    assert window.end == 20_000_000
    assert window.bed_start0() == 9_999_999
    assert window.bed_end() == 20_000_000


def test_intersect_bed_with_region(tmp_path: Path) -> None:
    bed = tmp_path / "truth.bed"
    bed.write_text(
        "\n".join(
            [
                "chr20\t0\t5000000",
                "chr20\t9000000\t15000000",
                "chr20\t19000000\t25000000",
                "chr1\t100\t200",
                "",
            ]
        ),
        encoding="utf-8",
    )
    out = tmp_path / "eval.bed"
    dry.intersect_bed_with_region(bed, "chr20:10000000-20000000", out)
    rows = [line.split("\t") for line in out.read_text(encoding="utf-8").splitlines() if line]
    # Region chr20:10000000-20000000 is BED [9999999, 20000000).
    assert rows == [
        ["chr20", "9999999", "15000000"],
        ["chr20", "19000000", "20000000"],
    ]


def test_downsample_view_command_construction() -> None:
    argv = dry.downsample_view_command(
        dry.ILLUMINA_BAM_URL,
        dry.REGION,
        dry.DOWNSAMPLE,
        dry.THREADS,
    )
    assert argv[:3] == ["samtools", "view", "-b"]
    assert "-s" in argv and argv[argv.index("-s") + 1] == "0.1"
    assert argv[-1] == dry.REGION
    assert argv[-2] == dry.ILLUMINA_BAM_URL


def test_minimap2_and_call_command_construction(tmp_path: Path) -> None:
    ref = tmp_path / "chr20.fa"
    index = tmp_path / "chr20.fa.mmi"
    r1 = tmp_path / "RO-abc_R1.fastq.gz"
    r2 = tmp_path / "RO-abc_R2.fastq.gz"
    ref.write_text(">chr20\nACGT\n", encoding="utf-8")
    index.write_text("idx", encoding="utf-8")
    r1.write_text("@r\nA\n+\nI\n", encoding="utf-8")
    r2.write_text("@r\nT\n+\nI\n", encoding="utf-8")

    align = dry.minimap2_align_command(ref, index, r1, r2, "RO-abc", 8)
    assert align[0] == "minimap2"
    assert "-ax" in align and "sr" in align
    assert any(arg.startswith("@RG") and "SM:RO-abc" in arg.replace("\\t", "\t") for arg in align)

    stages = dry.call_pipeline_commands(ref, dry.REGION, 4)
    assert stages[0][0] == "bcftools" and stages[0][1] == "mpileup"
    assert "-a" in stages[0] and stages[0][stages[0].index("-a") + 1] == "AD,DP"
    assert stages[1][1:] == ["call", "-mv", "-Ou", "--threads", "4"]
    assert stages[2][1] == "norm" and "-m" in stages[2] and "-any" in stages[2]


def test_vcfeval_command_construction(tmp_path: Path) -> None:
    sdf = tmp_path / "chr20.sdf"
    base = tmp_path / "truth.vcf.gz"
    calls = tmp_path / "calls.vcf.gz"
    bed = tmp_path / "eval.bed"
    out = tmp_path / "vcfeval_snv"
    argv = dry.vcfeval_command(sdf, base, calls, bed, out)
    assert argv[:2] == ["rtg", "vcfeval"]
    assert "-e" in argv and argv[argv.index("-e") + 1] == str(bed)
    assert "-o" in argv and argv[argv.index("-o") + 1] == str(out)


def test_parse_vcfeval_summary_fixture() -> None:
    text = _SUMMARY_FIXTURE.read_text(encoding="utf-8")
    summary = dry.parse_vcfeval_summary(text)
    assert summary.true_pos_baseline == 100
    assert summary.true_pos_call == 98
    assert summary.false_pos == 5
    assert summary.false_neg == 7
    assert summary.precision == pytest.approx(0.9515)
    assert summary.recall == pytest.approx(0.9346)
    assert summary.f1 == pytest.approx(0.9430)
