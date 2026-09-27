"""Intake report and custodian linkage table.

The report is not allowed to contain an original sample id. The linkage
table is the only file that pairs original ids with pseudonyms, and it is
written outside the repository.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from collections.abc import Sequence
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import TypedDict

from rogen_aging.intake.errors import IntakeError
from rogen_aging.intake.manifest import ManifestRow
from rogen_aging.intake.rewrite import ProducedFile


class _FileRecord(TypedDict):
    pseudonym: str
    file_type: str
    output_name: str
    checksum_type: str
    checksum: str
    pod5_metadata_scan: str | None


class IntakeReport(TypedDict):
    """JSON report written under the processing output directory."""

    timestamp: str
    git_commit: str
    leak_check: str
    checksums_verified: bool
    tool_versions: dict[str, str]
    files: list[_FileRecord]


def build_report(
    files: Sequence[ProducedFile],
    repo_root: Path,
) -> IntakeReport:
    """Build the intake report. Original ids are not accepted as arguments."""
    records: list[_FileRecord] = []
    for produced in files:
        records.append(
            {
                "pseudonym": produced.pseudonym,
                "file_type": produced.file_type,
                "output_name": produced.output_path.name,
                "checksum_type": produced.checksum_type,
                "checksum": produced.checksum,
                "pod5_metadata_scan": produced.pod5_metadata_scan,
            }
        )
    records.sort(key=lambda item: item["output_name"])
    return {
        "timestamp": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "git_commit": _git_commit(repo_root),
        "leak_check": "pass",
        "checksums_verified": True,
        "tool_versions": _tool_versions(),
        "files": records,
    }


def report_text(report: IntakeReport) -> str:
    """Serialize a report to JSON."""
    return json.dumps(report, indent=2, sort_keys=True) + "\n"


def assert_report_excludes_original_ids(text: str, original_ids: Sequence[str]) -> None:
    """Refuse to publish a report that still contains an original sample id.

    This check stays in force when Python is run with ``-O``.

    Args:
        text: Serialized report.
        original_ids: Original sample ids from the manifest.

    Raises:
        IntakeError: If any original id is a case-insensitive substring.
    """
    lowered = text.lower()
    for original_id in original_ids:
        if original_id.lower() in lowered:
            raise IntakeError("intake report contains an original sample id")
    assert not any(
        original_id.lower() in lowered for original_id in original_ids
    ), "intake report contains an original sample id"


def render_linkage(rows: Sequence[ManifestRow], pseudonyms: dict[str, str]) -> str:
    """Render the custodian linkage table."""
    seen: set[str] = set()
    lines = ["original_sample_id\tpseudonym"]
    for row in rows:
        if row.original_sample_id in seen:
            continue
        seen.add(row.original_sample_id)
        lines.append(f"{row.original_sample_id}\t{pseudonyms[row.original_sample_id]}")
    return "\n".join(lines) + "\n"


def _git_commit(repo_root: Path) -> str:
    proc = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
        timeout=10,
    )
    if proc.returncode != 0:
        return "unknown"
    commit = proc.stdout.strip()
    if not commit:
        return "unknown"
    return commit


def _tool_versions() -> dict[str, str]:
    return {
        "python": sys.version.split()[0],
        "samtools": _cli_version("samtools"),
        "bcftools": _cli_version("bcftools"),
        "pysam": _package_version("pysam"),
        "pod5": _package_version("pod5"),
    }


def _cli_version(command: str) -> str:
    if shutil.which(command) is None:
        return "not installed"
    proc = subprocess.run(
        [command, "--version"],
        capture_output=True,
        text=True,
        check=False,
        timeout=10,
    )
    text = (proc.stdout or proc.stderr).strip()
    if not text:
        return "unknown"
    return text.splitlines()[0].strip()


def _package_version(name: str) -> str:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return "not installed"
