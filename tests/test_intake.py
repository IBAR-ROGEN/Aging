"""Synthetic-fixture tests for sequencing intake.

Fixtures are created inside ``tmp_path``. No sequencing data is downloaded.
"""

from __future__ import annotations

import hashlib
import shutil
import stat
import subprocess
from pathlib import Path

import pytest
from typer.testing import CliRunner

from rogen_aging.cli.intake import app
from rogen_aging.config import find_repo_root
from rogen_aging.intake import IntakeError, make_key, pseudonym_for, run_intake
from rogen_aging.intake.key import KEY_ENV_VAR
from rogen_aging.intake.rewrite import rewrite_sam_header

_SAMPLE = "SMP-ALPHA-7Q"
_COLUMNS = (
    "original_sample_id",
    "file_name",
    "file_type",
    "checksum_type",
    "checksum",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    columns = list(_COLUMNS)
    if any("aliases" in row for row in rows):
        columns.append("aliases")
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(row.get(column, "") for column in columns))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _row(file_name: str, file_type: str, checksum: str) -> dict[str, str]:
    return {
        "original_sample_id": _SAMPLE,
        "file_name": file_name,
        "file_type": file_type,
        "checksum_type": "sha256",
        "checksum": checksum,
    }


def _place_key(path: Path, monkeypatch: pytest.MonkeyPatch) -> bytes:
    make_key(path)
    monkeypatch.setenv(KEY_ENV_VAR, str(path))
    return path.read_bytes()


def test_checksum_mismatch_writes_nothing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A bad checksum stops before any output or linkage file is created."""
    delivery = tmp_path / "delivery"
    delivery.mkdir()
    payload = delivery / "reads.fastq"
    payload.write_text("@READ\nACGT\n+\nFFFF\n", encoding="utf-8")
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(manifest, [_row("reads.fastq", "fastq", "0" * 64)])
    out_dir = tmp_path / "out"
    custodian = tmp_path / "custodian"
    _place_key(tmp_path / "pseudo.key", monkeypatch)

    with pytest.raises(IntakeError, match="checksum mismatch for reads.fastq"):
        run_intake(manifest, delivery, out_dir, custodian)

    assert not out_dir.exists()
    assert not custodian.exists()


def test_key_inside_repo_and_group_readable_key_are_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The key must sit outside the repo and must not be group- or world-readable."""
    repo = find_repo_root()
    forbidden = repo / "intake-test-do-not-commit.key"
    try:
        with pytest.raises(IntakeError, match="repository"):
            make_key(forbidden)
        assert not forbidden.exists()
    finally:
        if forbidden.exists():
            forbidden.unlink()

    delivery = tmp_path / "delivery"
    delivery.mkdir()
    payload = delivery / "reads.fastq"
    payload.write_text("@READ\nACGT\n+\nFFFF\n", encoding="utf-8")
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(manifest, [_row("reads.fastq", "fastq", _sha256(payload))])
    out_dir = tmp_path / "out"
    custodian = tmp_path / "custodian"

    monkeypatch.setenv(KEY_ENV_VAR, str(repo / "pyproject.toml"))
    with pytest.raises(IntakeError, match="repository"):
        run_intake(manifest, delivery, out_dir, custodian)
    assert not out_dir.exists()

    world_readable = tmp_path / "world.key"
    world_readable.write_bytes(b"k" * 32)
    world_readable.chmod(0o644)
    assert stat.S_IMODE(world_readable.stat().st_mode) & (stat.S_IRGRP | stat.S_IROTH)
    monkeypatch.setenv(KEY_ENV_VAR, str(world_readable))
    with pytest.raises(IntakeError, match="readable by group or others"):
        run_intake(manifest, delivery, out_dir, custodian)
    assert not out_dir.exists()
    assert not (custodian / "linkage.tsv").exists()


def test_pseudonym_is_stable_for_one_key_and_changes_with_the_key() -> None:
    """Same key and id always match; a different key does not."""
    first = b"a" * 32
    second = b"b" * 32
    again = pseudonym_for(first, _SAMPLE)
    assert again == pseudonym_for(first, _SAMPLE)
    assert again != pseudonym_for(second, _SAMPLE)
    assert again.startswith("RO-")
    assert len(again) == 15


def test_fastq_read_header_with_original_id_fails_leak_check(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """FASTQ is not rewritten; an original id in a read header deletes the output."""
    delivery = tmp_path / "delivery"
    delivery.mkdir()
    payload = delivery / "reads.fastq"
    payload.write_text(f"@{_SAMPLE}\nACGT\n+\nFFFF\n", encoding="utf-8")
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(manifest, [_row("reads.fastq", "fastq", _sha256(payload))])
    out_dir = tmp_path / "out"
    custodian = tmp_path / "custodian"
    _place_key(tmp_path / "pseudo.key", monkeypatch)

    with pytest.raises(IntakeError, match="FASTQ read header"):
        run_intake(manifest, delivery, out_dir, custodian)

    assert not out_dir.exists()
    assert not (custodian / "linkage.tsv").exists()


def test_vcf_sample_header_is_rewritten_and_report_has_no_original_id(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The VCF sample column becomes the pseudonym, and the report omits the original id."""
    delivery = tmp_path / "delivery"
    delivery.mkdir()
    payload = delivery / "calls.vcf"
    payload.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##source=synthetic-fixture",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + _SAMPLE,
                "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(manifest, [_row("calls.vcf", "vcf", _sha256(payload))])
    out_dir = tmp_path / "out"
    custodian = tmp_path / "custodian"
    key = _place_key(tmp_path / "pseudo.key", monkeypatch)

    run_intake(manifest, delivery, out_dir, custodian)

    pseudonym = pseudonym_for(key, _SAMPLE)
    rewritten = (out_dir / f"{pseudonym}.vcf").read_text(encoding="utf-8")
    chrom_line = next(line for line in rewritten.splitlines() if line.startswith("#CHROM"))
    assert chrom_line.split("\t")[-1] == pseudonym
    assert _SAMPLE not in rewritten
    report = (out_dir / "intake_report.json").read_text(encoding="utf-8")
    assert _SAMPLE.lower() not in report.lower()
    linkage = (custodian / "linkage.tsv").read_text(encoding="utf-8")
    assert f"{_SAMPLE}\t{pseudonym}" in linkage


def test_bam_read_group_rewrite(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Rewrite @RG when samtools is available; otherwise skip."""
    if shutil.which("samtools") is None:
        pytest.skip("samtools is not installed")

    delivery = tmp_path / "delivery"
    delivery.mkdir()
    sam = delivery / "reads.sam"
    sam.write_text(
        "\n".join(
            [
                "@HD\tVN:1.6\tSO:unsorted",
                "@SQ\tSN:chr1\tLN:1000",
                f"@RG\tID:{_SAMPLE}\tSM:{_SAMPLE}\tLB:{_SAMPLE}",
                f"@PG\tID:bwa\tPN:bwa\tCL:bwa mem ref {_SAMPLE}.fastq",
                f"r1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tF\tRG:Z:{_SAMPLE}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    bam = delivery / "reads.bam"
    built = shutil.which("samtools")
    assert built is not None
    proc = subprocess.run(
        [built, "view", "-b", "-o", str(bam), str(sam)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(manifest, [_row("reads.bam", "bam", _sha256(bam))])
    out_dir = tmp_path / "out"
    custodian = tmp_path / "custodian"
    key = _place_key(tmp_path / "pseudo.key", monkeypatch)

    run_intake(manifest, delivery, out_dir, custodian)

    pseudonym = pseudonym_for(key, _SAMPLE)
    header = subprocess.run(
        [built, "view", "-H", str(out_dir / f"{pseudonym}.bam")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert header.returncode == 0, header.stderr
    assert _SAMPLE not in header.stdout
    assert f"SM:{pseudonym}" in header.stdout
    report = (out_dir / "intake_report.json").read_text(encoding="utf-8")
    assert _SAMPLE.lower() not in report.lower()


def test_pod5_without_package_is_not_described_as_clean(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A renamed POD5 file is not called clean when its metadata was not scanned."""
    import importlib.util

    if importlib.util.find_spec("pod5") is not None:
        pytest.skip("pod5 is installed; this fixture only covers the unscanned path")
    delivery = tmp_path / "delivery"
    delivery.mkdir()
    payload = delivery / "reads.pod5"
    payload.write_bytes(b"POD5-fixture-not-a-real-file")
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(manifest, [_row("reads.pod5", "pod5", _sha256(payload))])
    out_dir = tmp_path / "out"
    custodian = tmp_path / "custodian"
    _place_key(tmp_path / "pseudo.key", monkeypatch)

    run_intake(manifest, delivery, out_dir, custodian)

    report = (out_dir / "intake_report.json").read_text(encoding="utf-8")
    assert "pod5 metadata not scanned" in report
    assert "clean" not in report.lower()
    assert _SAMPLE.lower() not in report.lower()


def test_alias_in_fastq_header_fails_leak_check(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An alias that is not the primary id is still a leak."""
    delivery = tmp_path / "delivery"
    delivery.mkdir()
    payload = delivery / "reads.fastq"
    payload.write_text("@NA24385\nACGT\n+\nFFFF\n", encoding="utf-8")
    manifest = tmp_path / "manifest.tsv"
    row = _row("reads.fastq", "fastq", _sha256(payload))
    row["original_sample_id"] = "HG002"
    row["aliases"] = "NA24385,GM24385"
    _write_manifest(manifest, [row])
    out_dir = tmp_path / "out"
    custodian = tmp_path / "custodian"
    _place_key(tmp_path / "pseudo.key", monkeypatch)

    with pytest.raises(IntakeError, match="FASTQ read header"):
        run_intake(manifest, delivery, out_dir, custodian)

    assert not out_dir.exists()
    assert not (custodian / "linkage.tsv").exists()


def test_rg_sm_alias_is_rewritten_in_the_header() -> None:
    """@RG SM may be an alias. That string leaves the header, and so does a @PG path."""
    header = "\n".join(
        [
            "@HD\tVN:1.6\tSO:unsorted",
            "@SQ\tSN:chr20\tLN:64444167",
            "@RG\tID:HG002\tSM:NA24385\tLB:GM24385",
            "@PG\tID:novoalign\tPN:novoalign\tCL:novoalign -d ref NA24385_R1.fastq",
            "",
        ]
    )
    rewritten = rewrite_sam_header(header, ("HG002", "NA24385", "GM24385"), "RO-0123456789ab")
    assert "NA24385" not in rewritten
    assert "GM24385" not in rewritten
    assert "HG002" not in rewritten
    assert "SM:RO-0123456789ab" in rewritten
    assert "@PG" not in rewritten


def test_bam_read_group_sm_alias_is_rewritten(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Intake removes an alias stored in @RG SM when samtools can rewrite the BAM."""
    if shutil.which("samtools") is None:
        pytest.skip("samtools is not installed")

    delivery = tmp_path / "delivery"
    delivery.mkdir()
    sam = delivery / "reads.sam"
    sam.write_text(
        "\n".join(
            [
                "@HD\tVN:1.6\tSO:unsorted",
                "@SQ\tSN:chr20\tLN:64444167",
                "@RG\tID:HG002\tSM:NA24385\tLB:GM24385",
                "@PG\tID:novoalign\tPN:novoalign\tCL:novoalign -d ref NA24385.fastq",
                "r1\t0\tchr20\t1\t60\t1M\t*\t0\t0\tA\tF\tRG:Z:HG002",
                "",
            ]
        ),
        encoding="utf-8",
    )
    bam = delivery / "reads.bam"
    built = shutil.which("samtools")
    assert built is not None
    proc = subprocess.run(
        [built, "view", "-b", "-o", str(bam), str(sam)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    manifest = tmp_path / "manifest.tsv"
    row = _row("reads.bam", "bam", _sha256(bam))
    row["original_sample_id"] = "HG002"
    row["aliases"] = "NA24385,GM24385"
    _write_manifest(manifest, [row])
    out_dir = tmp_path / "out"
    custodian = tmp_path / "custodian"
    key = _place_key(tmp_path / "pseudo.key", monkeypatch)

    run_intake(manifest, delivery, out_dir, custodian)

    pseudonym = pseudonym_for(key, "HG002")
    header = subprocess.run(
        [built, "view", "-H", str(out_dir / f"{pseudonym}.bam")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert header.returncode == 0, header.stderr
    assert "NA24385" not in header.stdout
    assert "GM24385" not in header.stdout
    assert "HG002" not in header.stdout
    assert f"SM:{pseudonym}" in header.stdout


def test_make_key_cli_does_not_print_the_key(tmp_path: Path) -> None:
    """``rogen-intake make-key`` writes 32 bytes and does not echo them."""
    dest = tmp_path / "cohort.key"
    result = CliRunner().invoke(app, ["make-key", "--out", str(dest)])
    assert result.exit_code == 0, result.output
    key = dest.read_bytes()
    assert len(key) == 32
    assert stat.S_IMODE(dest.stat().st_mode) == 0o600
    assert key.hex() not in result.output
