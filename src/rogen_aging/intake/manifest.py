"""Delivery manifest parsing and checksum verification."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from rogen_aging.intake.errors import IntakeError

MANIFEST_COLUMNS: tuple[str, ...] = (
    "original_sample_id",
    "file_name",
    "file_type",
    "checksum_type",
    "checksum",
)
OPTIONAL_MANIFEST_COLUMNS: tuple[str, ...] = ("aliases",)
ALLOWED_FILE_TYPES: frozenset[str] = frozenset({"fastq", "bam", "cram", "vcf", "pod5"})
ALLOWED_CHECKSUM_TYPES: frozenset[str] = frozenset({"md5", "sha256"})


@dataclass(frozen=True)
class ManifestRow:
    """One file in a delivery manifest."""

    original_sample_id: str
    file_name: str
    file_type: str
    checksum_type: str
    checksum: str
    aliases: tuple[str, ...] = ()

    def search_ids(self) -> tuple[str, ...]:
        """Original id plus aliases, in order, without case-insensitive duplicates."""
        ordered: list[str] = []
        seen: set[str] = set()
        for sample_id in (self.original_sample_id, *self.aliases):
            key = sample_id.casefold()
            if not sample_id or key in seen:
                continue
            seen.add(key)
            ordered.append(sample_id)
        return tuple(ordered)


def load_manifest(path: Path) -> list[ManifestRow]:
    """Load a TSV manifest and reject unexpected columns or file types.

    Args:
        path: Path to the delivery manifest.

    Returns:
        One row per file, in file order.

    Raises:
        IntakeError: If the header, a file type, or a checksum type is invalid.
    """
    if not path.is_file():
        raise IntakeError(f"manifest not found: {path}")
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    columns = [str(name) for name in frame.columns.tolist()]
    allowed = [
        list(MANIFEST_COLUMNS),
        [*MANIFEST_COLUMNS, *OPTIONAL_MANIFEST_COLUMNS],
    ]
    if columns not in allowed:
        raise IntakeError(
            "manifest columns must be "
            + ", ".join(MANIFEST_COLUMNS)
            + " with optional aliases; found "
            + ", ".join(columns)
        )
    rows: list[ManifestRow] = []
    records = frame.to_dict(orient="records")
    for index, record in enumerate(records, start=1):
        values = {str(key): str(value).strip() for key, value in record.items()}
        sample_id = values["original_sample_id"]
        file_name = values["file_name"]
        file_type = values["file_type"].lower()
        checksum_type = values["checksum_type"].lower()
        checksum = values["checksum"].lower()
        if not sample_id or not file_name or not checksum:
            raise IntakeError(f"manifest row {index} has an empty required field")
        if file_type not in ALLOWED_FILE_TYPES:
            raise IntakeError(
                f"manifest file_type {file_type!r} is outside fastq|bam|cram|vcf|pod5"
            )
        if checksum_type not in ALLOWED_CHECKSUM_TYPES:
            raise IntakeError(f"manifest checksum_type {checksum_type!r} is outside md5|sha256")
        rows.append(
            ManifestRow(
                original_sample_id=sample_id,
                file_name=file_name,
                file_type=file_type,
                checksum_type=checksum_type,
                checksum=checksum,
                aliases=_parse_aliases(values.get("aliases", ""), sample_id),
            )
        )
    if not rows:
        raise IntakeError("manifest has no file rows")
    return rows


def _parse_aliases(raw: str, primary: str) -> tuple[str, ...]:
    """Split a comma-separated alias cell. The primary id is not repeated."""
    found: list[str] = []
    seen = {primary.casefold()}
    for part in raw.split(","):
        alias = part.strip()
        key = alias.casefold()
        if not alias or key in seen:
            continue
        seen.add(key)
        found.append(alias)
    return tuple(found)


def identifiers_to_scrub(rows: list[ManifestRow]) -> list[str]:
    """Every original id and alias the leak check must search."""
    found: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for sample_id in row.search_ids():
            key = sample_id.casefold()
            if key in seen:
                continue
            seen.add(key)
            found.append(sample_id)
    return found


def delivery_file(delivery_dir: Path, file_name: str) -> Path:
    """Resolve ``file_name`` inside ``delivery_dir``.

    Args:
        delivery_dir: Directory that holds the raw delivery.
        file_name: Manifest file name, relative to ``delivery_dir``.

    Returns:
        Absolute path of the delivery file.

    Raises:
        IntakeError: If the name is absolute or escapes ``delivery_dir``.
    """
    raw = Path(file_name)
    if raw.is_absolute() or ".." in raw.parts:
        raise IntakeError(f"refusing delivery path outside the delivery directory: {file_name}")
    root = delivery_dir.resolve()
    candidate = (root / raw).resolve()
    if not candidate.is_relative_to(root):
        raise IntakeError(f"refusing delivery path outside the delivery directory: {file_name}")
    return candidate


def _digest(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1 << 16)
            if not chunk:
                break
            hasher.update(chunk)
    return hasher.hexdigest()


def verify_checksums(delivery_dir: Path, rows: list[ManifestRow]) -> None:
    """Check every manifest checksum. Write nothing.

    Args:
        delivery_dir: Directory that holds the raw delivery.
        rows: Parsed manifest rows.

    Raises:
        IntakeError: On the first missing file or checksum mismatch.
            The message names that file.
    """
    if not delivery_dir.is_dir():
        raise IntakeError(f"delivery directory not found: {delivery_dir}")
    for row in rows:
        path = delivery_file(delivery_dir, row.file_name)
        if not path.is_file():
            raise IntakeError(f"missing file: {row.file_name}")
        actual = _digest(path, row.checksum_type)
        if actual.lower() != row.checksum.lower():
            raise IntakeError(f"checksum mismatch for {row.file_name}")
