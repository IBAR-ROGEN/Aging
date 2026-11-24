"""Synthetic BAM tests for relative mtDNA copy-number estimation.

BAMs are built inside this module and labelled synthetic. Expected estimates
are computed by hand from the planted depths, not from literature values.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pysam
import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SCRIPT = _REPO_ROOT / "scripts" / "mtdna" / "estimate_mtdna_copy_number.py"


def _load_script() -> ModuleType:
    name = "estimate_mtdna_copy_number"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


mod = _load_script()

# Contig length (bp) for every synthetic chromosome, including chrM.
_SYNTHETIC_CONTIG_LEN = 50
# Full-length reads per autosome → meandepth = 2 on each of chr1–chr22.
_AUTO_READS = 2
# Full-length reads on chrM → meandepth = 10.
_MT_READS = 10
# Hand calculation:
#   mean_chrM = 10
#   mean_chr1_22 = 2  (all autosomes equal length and depth)
#   estimate = 2 * (10 / 2) = 10
_HAND_EXPECTED_ESTIMATE = 10.0


def _write_synthetic_bam(
    path: Path,
    *,
    auto_reads: int,
    mt_reads: int,
    contig_len: int = _SYNTHETIC_CONTIG_LEN,
) -> Path:
    """Write a tiny coordinate-sorted synthetic BAM with chr1–chr22 and chrM."""
    contigs = [f"chr{i}" for i in range(1, 23)] + ["chrM"]
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": contig_len} for name in contigs],
        "RG": [{"ID": "synthetic", "SM": "synthetic", "DS": "synthetic fixture"}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for chrom_idx, chrom in enumerate(contigs):
            n_reads = mt_reads if chrom == "chrM" else auto_reads
            for read_idx in range(n_reads):
                aligned = pysam.AlignedSegment()
                aligned.query_name = f"synthetic_{chrom}_{read_idx}"
                aligned.query_sequence = "A" * contig_len
                aligned.flag = 0
                aligned.reference_id = chrom_idx
                aligned.reference_start = 0
                aligned.mapping_quality = 60
                aligned.cigar = ((0, contig_len),)
                aligned.query_qualities = pysam.qualitystring_to_array("I" * contig_len)
                aligned.set_tag("RG", "synthetic")
                out.write(aligned)
    pysam.index(str(path))
    return path


@pytest.fixture(scope="module")
def samtools_path() -> str:
    """Resolve samtools or fail loudly (do not substitute another method)."""
    return mod.find_samtools()


def test_warning_constant_matches_required_text() -> None:
    assert mod.WARNING == (
        "Blood-derived estimates are confounded by platelets and cell composition; "
        "values are relative, not absolute."
    )


def test_synthetic_bam_hand_computed_estimate(tmp_path: Path, samtools_path: str) -> None:
    """Planted depths: chrM=10, each autosome=2 → estimate = 2*(10/2) = 10."""
    bam = _write_synthetic_bam(
        tmp_path / "synthetic_mtdna_cn.bam",
        auto_reads=_AUTO_READS,
        mt_reads=_MT_READS,
    )
    metrics = mod.estimate_one_bam(samtools_path, bam)
    assert metrics["mean_chrM_depth"] == pytest.approx(float(_MT_READS))
    assert metrics["mean_depth_chr1_22"] == pytest.approx(float(_AUTO_READS))
    assert metrics["mtdna_copy_number_estimate"] == pytest.approx(_HAND_EXPECTED_ESTIMATE)


def test_output_header_states_platelet_warning(tmp_path: Path, samtools_path: str) -> None:
    bam = _write_synthetic_bam(
        tmp_path / "synthetic_mtdna_cn_header.bam",
        auto_reads=_AUTO_READS,
        mt_reads=_MT_READS,
    )
    samples_csv = tmp_path / "synthetic_samples.csv"
    samples_csv.write_text(
        "sample_id,bam_path\n" f"synthetic_A,{bam.name}\n",
        encoding="utf-8",
    )

    out_csv = tmp_path / "out" / "estimates.csv"
    out_prov = tmp_path / "out" / "estimates_provenance.json"
    frame = mod.run_estimates(samples_csv, out_csv, out_prov)

    text = out_csv.read_text(encoding="utf-8")
    first_line = text.splitlines()[0]
    assert first_line == f"# {mod.WARNING}"
    assert "synthetic_A" in text
    assert frame.height == 1
    assert frame["mtdna_copy_number_estimate"][0] == pytest.approx(_HAND_EXPECTED_ESTIMATE)

    provenance = out_prov.read_text(encoding="utf-8")
    assert mod.WARNING in provenance
    assert "synthetic" in provenance.lower()
    assert samtools_path in provenance


def test_missing_samples_csv_fails_loudly(tmp_path: Path) -> None:
    missing = tmp_path / "does_not_exist.csv"
    with pytest.raises(AssertionError, match="MISSING INPUT"):
        mod.load_samples(missing)


def test_second_synthetic_ratio(tmp_path: Path, samtools_path: str) -> None:
    """chrM=6, autosomes=3 → estimate = 2*(6/3) = 4 (hand-computed)."""
    bam = _write_synthetic_bam(
        tmp_path / "synthetic_mtdna_cn_ratio2.bam",
        auto_reads=3,
        mt_reads=6,
    )
    metrics = mod.estimate_one_bam(samtools_path, bam)
    assert metrics["mtdna_copy_number_estimate"] == pytest.approx(4.0)
