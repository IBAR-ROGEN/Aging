"""Rewrite delivery files so output names and headers use pseudonyms.

FASTQ reads are copied unchanged. POD5 files are copied unchanged.
BAM/CRAM headers go through ``samtools reheader`` when samtools is on
``PATH``, otherwise through pysam when that package is installed.
VCF sample columns go through ``bcftools reheader --samples`` when
bcftools is on ``PATH``. A text rewrite of the sample column is used
only when bcftools is absent, so a fixture VCF can still be checked.
"""

from __future__ import annotations

import gzip
import importlib
import importlib.util
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rogen_aging.intake.errors import IntakeError
from rogen_aging.intake.manifest import ManifestRow, delivery_file

_COMPOUND_SUFFIXES: tuple[str, ...] = (".fastq.gz", ".fq.gz", ".vcf.gz", ".vcf.bgz")
_RG_FIELDS = frozenset({"ID", "SM", "LB"})


@dataclass
class ProducedFile:
    """One pseudonymized output file, before it is moved into place."""

    pseudonym: str
    file_type: str
    output_path: Path
    checksum_type: str
    checksum: str
    pod5_metadata_scan: str | None = None


def output_filename(pseudonym: str, file_name: str, file_type: str) -> str:
    """Return ``<pseudonym>.<original extension>``."""
    name = Path(file_name).name
    lower = name.lower()
    for suffix in _COMPOUND_SUFFIXES:
        if lower.endswith(suffix):
            return f"{pseudonym}{name[-len(suffix) :]}"
    suffix = Path(name).suffix
    if not suffix:
        suffix = f".{file_type}"
    return f"{pseudonym}{suffix}"


def rewrite_delivery(
    work_dir: Path,
    delivery_dir: Path,
    rows: list[ManifestRow],
    pseudonyms: dict[str, str],
) -> list[ProducedFile]:
    """Write pseudonymized copies under ``work_dir``.

    Args:
        work_dir: Temporary directory. Not the final output directory.
        delivery_dir: Raw delivery directory.
        rows: Manifest rows whose checksums have already been verified.
        pseudonyms: Map of original sample id to pseudonym.

    Returns:
        Produced files inside ``work_dir``.

    Raises:
        IntakeError: If a required external tool is missing or a rewrite fails.
    """
    produced: list[ProducedFile] = []
    used_names: set[str] = set()
    for row in rows:
        pseudonym = pseudonyms[row.original_sample_id]
        name = output_filename(pseudonym, row.file_name, row.file_type)
        if name in used_names:
            raise IntakeError(f"output file name collision for {name}")
        used_names.add(name)
        source = delivery_file(delivery_dir, row.file_name)
        target = work_dir / name
        target.parent.mkdir(parents=True, exist_ok=True)
        if row.file_type in {"fastq", "pod5"}:
            shutil.copyfile(source, target)
        elif row.file_type in {"bam", "cram"}:
            _rewrite_alignment(source, target, row.search_ids(), pseudonym, work_dir)
        elif row.file_type == "vcf":
            _rewrite_vcf(source, target, row.search_ids(), pseudonym, work_dir)
        else:
            raise IntakeError(f"unsupported file_type {row.file_type!r}")
        produced.append(
            ProducedFile(
                pseudonym=pseudonym,
                file_type=row.file_type,
                output_path=target,
                checksum_type=row.checksum_type,
                checksum=row.checksum,
            )
        )
    return produced


def _rewrite_alignment(
    source: Path,
    target: Path,
    original_ids: tuple[str, ...],
    pseudonym: str,
    work_dir: Path,
) -> None:
    if shutil.which("samtools"):
        _samtools_reheader(source, target, original_ids, pseudonym, work_dir)
        return
    pysam_mod = _optional_module("pysam")
    if pysam_mod is None:
        raise IntakeError(
            "samtools is not installed and pysam is not installed; cannot rewrite BAM or CRAM"
        )
    _pysam_reheader(pysam_mod, source, target, original_ids, pseudonym)


def _samtools_reheader(
    source: Path,
    target: Path,
    original_ids: tuple[str, ...],
    pseudonym: str,
    work_dir: Path,
) -> None:
    viewed = subprocess.run(
        ["samtools", "view", "-H", str(source)],
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )
    if viewed.returncode != 0:
        detail = (viewed.stderr or "").strip()
        raise IntakeError(f"samtools view -H failed: {detail}")
    header_text = rewrite_sam_header(viewed.stdout, original_ids, pseudonym)
    header_path = work_dir / f".header-{pseudonym}.sam"
    header_path.write_text(header_text, encoding="utf-8")
    with target.open("wb") as handle:
        rewritten = subprocess.run(
            ["samtools", "reheader", str(header_path), str(source)],
            stdout=handle,
            stderr=subprocess.PIPE,
            check=False,
            timeout=60,
        )
    header_path.unlink(missing_ok=True)
    if rewritten.returncode != 0:
        detail = rewritten.stderr.decode("utf-8", errors="replace").strip()
        raise IntakeError(f"samtools reheader failed: {detail}")
    _rewrite_rg_tags(target, original_ids, pseudonym, work_dir)


def _rewrite_rg_tags(
    bam_or_cram: Path,
    original_ids: tuple[str, ...],
    pseudonym: str,
    work_dir: Path,
) -> None:
    """Replace alignment ``RG:Z`` values that still carry the original id.

    ``samtools reheader`` updates the header only. Read-group tags on
    alignments are rewritten here so the original id is not left in the file.
    """
    streamed = subprocess.run(
        ["samtools", "view", "-h", str(bam_or_cram)],
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )
    if streamed.returncode != 0:
        detail = (streamed.stderr or "").strip()
        raise IntakeError(f"samtools view -h failed: {detail}")
    lines: list[str] = []
    replacement = f"RG:Z:{pseudonym}"
    needles = {f"rg:z:{sample_id.casefold()}" for sample_id in original_ids}
    for line in streamed.stdout.splitlines():
        if line.startswith("@") or "RG:Z:" not in line.upper():
            lines.append(line)
            continue
        parts = []
        for field in line.split("\t"):
            if field.casefold() in needles:
                parts.append(replacement)
            else:
                parts.append(field)
        lines.append("\t".join(parts))
    sam_text = "\n".join(lines) + "\n"
    sam_path = work_dir / f".body-{pseudonym}.sam"
    sam_path.write_text(sam_text, encoding="utf-8")
    rewritten = bam_or_cram.with_suffix(bam_or_cram.suffix + ".rgfix")
    proc = subprocess.run(
        ["samtools", "view", "-b", "-o", str(rewritten), str(sam_path)],
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )
    sam_path.unlink(missing_ok=True)
    if proc.returncode != 0:
        detail = (proc.stderr or "").strip()
        raise IntakeError(f"samtools view -b failed while rewriting RG tags: {detail}")
    rewritten.replace(bam_or_cram)


def rewrite_sam_header(
    header: str, original_ids: tuple[str, ...] | list[str], pseudonym: str
) -> str:
    """Rewrite ``@RG`` ID/SM/LB and drop ``@PG`` lines whose CL contains any id."""
    lines: list[str] = []
    for raw_line in header.splitlines():
        if raw_line.startswith("@PG") and _pg_cl_contains(raw_line, original_ids):
            continue
        if raw_line.startswith("@RG"):
            lines.append(_rewrite_rg_line(raw_line, pseudonym))
            continue
        lines.append(_replace_ids(raw_line, original_ids, pseudonym))
    text = "\n".join(lines)
    if text and not text.endswith("\n"):
        text += "\n"
    return text


def _pg_cl_contains(line: str, original_ids: tuple[str, ...] | list[str]) -> bool:
    for field in line.split("\t"):
        tag, sep, value = field.partition(":")
        if tag == "CL" and sep and _contains_id(value, original_ids):
            return True
    return False


def _rewrite_rg_line(line: str, pseudonym: str) -> str:
    parts = line.split("\t")
    rewritten = [parts[0]]
    seen: set[str] = set()
    for field in parts[1:]:
        tag, sep, _value = field.partition(":")
        if sep and tag in _RG_FIELDS:
            rewritten.append(f"{tag}:{pseudonym}")
            seen.add(tag)
        else:
            rewritten.append(field)
    for tag in ("ID", "SM", "LB"):
        if tag not in seen:
            rewritten.append(f"{tag}:{pseudonym}")
    return "\t".join(rewritten)


def _replace_ids(text: str, original_ids: tuple[str, ...] | list[str], pseudonym: str) -> str:
    for original_id in sorted(original_ids, key=len, reverse=True):
        text = re.sub(re.escape(original_id), pseudonym, text, flags=re.IGNORECASE)
    return text


def _contains_id(text: str, original_ids: tuple[str, ...] | list[str]) -> bool:
    lowered = text.lower()
    return any(sample_id.lower() in lowered for sample_id in original_ids)


def _pysam_reheader(
    pysam_mod: Any,
    source: Path,
    target: Path,
    original_ids: tuple[str, ...],
    pseudonym: str,
) -> None:
    read_mode = "rc" if source.suffix.lower() == ".cram" else "rb"
    write_mode = "wc" if target.suffix.lower() == ".cram" else "wb"
    with pysam_mod.AlignmentFile(str(source), read_mode) as inbound:
        header = inbound.header.to_dict()
        for read_group in header.get("RG", []):
            read_group["ID"] = pseudonym
            read_group["SM"] = pseudonym
            read_group["LB"] = pseudonym
        kept_programs = []
        for program in header.get("PG", []):
            command = str(program.get("CL", ""))
            if _contains_id(command, original_ids):
                continue
            kept_programs.append(program)
        header["PG"] = kept_programs
        with pysam_mod.AlignmentFile(str(target), write_mode, header=header) as outbound:
            for record in inbound:
                tag = record.get_tag("RG") if record.has_tag("RG") else None
                if isinstance(tag, str) and _contains_id(tag, original_ids):
                    record.set_tag("RG", pseudonym, value_type="Z")
                outbound.write(record)


def _rewrite_vcf(
    source: Path,
    target: Path,
    original_ids: tuple[str, ...],
    pseudonym: str,
    work_dir: Path,
) -> None:
    if shutil.which("bcftools"):
        _bcftools_reheader(source, target, pseudonym, work_dir)
        _scrub_vcf_header_id(target, original_ids, pseudonym)
        return
    _python_vcf_reheader(source, target, original_ids, pseudonym)


def _bcftools_reheader(source: Path, target: Path, pseudonym: str, work_dir: Path) -> None:
    samples_path = work_dir / f".samples-{pseudonym}.txt"
    sample_count = _vcf_sample_count(source)
    samples_path.write_text("\n".join([pseudonym] * sample_count) + "\n", encoding="utf-8")
    proc = subprocess.run(
        [
            "bcftools",
            "reheader",
            "--samples",
            str(samples_path),
            "-o",
            str(target),
            str(source),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )
    samples_path.unlink(missing_ok=True)
    if proc.returncode != 0:
        detail = (proc.stderr or "").strip()
        raise IntakeError(f"bcftools reheader failed: {detail}")


def _scrub_vcf_header_id(path: Path, original_ids: tuple[str, ...], pseudonym: str) -> None:
    """Replace any original id still present in VCF header lines after bcftools."""
    compressed = _is_gzip(path)
    text = _read_text(path)
    lines = text.splitlines(keepends=True)
    changed = False
    rewritten: list[str] = []
    for line in lines:
        if line.startswith("#") and _contains_id(line, original_ids):
            rewritten.append(_replace_ids(line, original_ids, pseudonym))
            changed = True
        else:
            rewritten.append(line)
    if changed:
        _write_text(path, "".join(rewritten), compressed)


def _python_vcf_reheader(
    source: Path, target: Path, original_ids: tuple[str, ...], pseudonym: str
) -> None:
    compressed = _is_gzip(source)
    text = _read_text(source)
    if "\0" in text:
        raise IntakeError("bcftools is not installed; cannot rewrite a binary VCF sample header")
    lines = text.splitlines(keepends=True)
    found_chrom = False
    rewritten: list[str] = []
    for line in lines:
        if line.startswith("#CHROM"):
            found_chrom = True
            newline = "\n" if line.endswith("\n") else ""
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                raise IntakeError("VCF header has no sample column")
            for index in range(9, len(parts)):
                parts[index] = pseudonym
            line = "\t".join(parts) + newline
        elif line.startswith("#"):
            line = _replace_ids(line, original_ids, pseudonym)
        rewritten.append(line)
    if not found_chrom:
        raise IntakeError("VCF is missing a #CHROM header line")
    _write_text(target, "".join(rewritten), compressed)


def _vcf_sample_count(path: Path) -> int:
    text = _read_text(path)
    for line in text.splitlines():
        if line.startswith("#CHROM"):
            parts = line.split("\t")
            if len(parts) < 10:
                raise IntakeError("VCF header has no sample column")
            return len(parts) - 9
    raise IntakeError("VCF is missing a #CHROM header line")


def _read_text(path: Path) -> str:
    if _is_gzip(path):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return handle.read()
    return path.read_text(encoding="utf-8", errors="replace")


def _write_text(path: Path, text: str, compressed: bool) -> None:
    if compressed:
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
        return
    path.write_text(text, encoding="utf-8")


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def _optional_module(name: str) -> Any | None:
    if importlib.util.find_spec(name) is None:
        return None
    return importlib.import_module(name)
