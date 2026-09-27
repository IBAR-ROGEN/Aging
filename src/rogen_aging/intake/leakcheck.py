"""Search pseudonymized outputs for original sample identifiers."""

from __future__ import annotations

import gzip
import importlib
import importlib.util
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from rogen_aging.intake.errors import IntakeError
from rogen_aging.intake.rewrite import ProducedFile

FASTQ_HEADER_LIMIT = 10_000
POD5_NOT_SCANNED = "pod5 metadata not scanned"


def leak_check(files: Sequence[ProducedFile], original_ids: Sequence[str]) -> None:
    """Fail closed when any original id remains in the searched surfaces.

    Searches output file names, BAM/CRAM/VCF headers, and the first 10,000
    FASTQ read headers. POD5 run metadata is searched only when the ``pod5``
    package is installed.

    Args:
        files: Pseudonymized files still in the temporary directory.
        original_ids: Every original sample id from the manifest.

    Raises:
        IntakeError: On the first hit. The message does not repeat the id.
    """
    ids = [sample_id for sample_id in original_ids if sample_id]
    for produced in files:
        if _contains_any(produced.output_path.name, ids):
            raise IntakeError(
                "leak check failed: an output file name contains an original sample id"
            )
        if produced.file_type in {"bam", "cram"}:
            header = _alignment_header(produced.output_path)
            if _contains_any(header, ids):
                raise IntakeError(
                    "leak check failed for "
                    f"{produced.output_path.name}: original sample id remains in the file header"
                )
        elif produced.file_type == "vcf":
            header = _vcf_header(produced.output_path)
            if _contains_any(header, ids):
                raise IntakeError(
                    "leak check failed for "
                    f"{produced.output_path.name}: original sample id remains in the VCF header"
                )
        elif produced.file_type == "fastq":
            for read_header in _fastq_headers(produced.output_path, FASTQ_HEADER_LIMIT):
                if _contains_any(read_header, ids):
                    raise IntakeError(
                        "leak check failed for "
                        f"{produced.output_path.name}: original sample id remains in a FASTQ read header"
                    )
        elif produced.file_type == "pod5":
            produced.pod5_metadata_scan = _scan_pod5(produced.output_path, ids)


def _contains_any(haystack: str, original_ids: Sequence[str]) -> bool:
    lowered = haystack.lower()
    return any(sample_id.lower() in lowered for sample_id in original_ids)


def _alignment_header(path: Path) -> str:
    if shutil.which("samtools"):
        proc = subprocess.run(
            ["samtools", "view", "-H", str(path)],
            capture_output=True,
            text=True,
            check=False,
            timeout=60,
        )
        if proc.returncode != 0:
            detail = (proc.stderr or "").strip()
            raise IntakeError(f"could not read alignment header for leak check: {detail}")
        return proc.stdout
    pysam_mod = _optional_module("pysam")
    if pysam_mod is None:
        raise IntakeError("could not read alignment header for leak check")
    mode = "rc" if path.suffix.lower() == ".cram" else "rb"
    with pysam_mod.AlignmentFile(str(path), mode) as handle:
        text = str(handle.header)
    return text


def _vcf_header(path: Path) -> str:
    lines: list[str] = []
    for line in _text_lines(path):
        if not line.startswith("#"):
            break
        lines.append(line)
    return "".join(lines)


def _fastq_headers(path: Path, limit: int) -> list[str]:
    headers: list[str] = []
    lines = _text_lines(path)
    index = 0
    while index < len(lines) and len(headers) < limit:
        header = lines[index]
        headers.append(header.rstrip("\n"))
        index += 4
    return headers


def _text_lines(path: Path) -> list[str]:
    if _is_gzip(path):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return handle.readlines()
    return path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def _scan_pod5(path: Path, original_ids: Sequence[str]) -> str:
    if importlib.util.find_spec("pod5") is None:
        return POD5_NOT_SCANNED
    module = importlib.import_module("pod5")
    reader_cls = getattr(module, "Reader", None)
    if reader_cls is None:
        return POD5_NOT_SCANNED
    texts: list[str] = []
    with reader_cls(str(path)) as reader:
        run_infos = getattr(reader, "run_infos", None)
        if run_infos is None:
            return POD5_NOT_SCANNED
        for info in run_infos:
            texts.extend(_stringify(info))
    if _contains_any("\n".join(texts), original_ids):
        raise IntakeError(
            "leak check failed for " f"{path.name}: original sample id remains in POD5 run metadata"
        )
    return "scanned"


def _stringify(value: Any, depth: int = 0) -> list[str]:
    if depth > 6:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, dict):
        found: list[str] = []
        for key, item in value.items():
            found.extend(_stringify(key, depth + 1))
            found.extend(_stringify(item, depth + 1))
        return found
    if isinstance(value, (list, tuple)):
        found = []
        for item in value:
            found.extend(_stringify(item, depth + 1))
        return found
    if hasattr(value, "__dict__"):
        return _stringify(vars(value), depth + 1)
    return []


def _optional_module(name: str) -> Any | None:
    if importlib.util.find_spec(name) is None:
        return None
    return importlib.import_module(name)
