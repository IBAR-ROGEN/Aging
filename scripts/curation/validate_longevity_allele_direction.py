#!/usr/bin/env python3
"""Validate ``longevity_allele_direction.csv`` completeness for the 47-panel.

Reports how many of 47 rows are complete, lists incomplete rsIDs, and flags
A/T and C/G palindromes. Allowed categorical values are taken only from
``docs/PENDING_ALLELES.md`` or documentation next to the curation file. If no
allow-list is documented, the script stops and lists distinct values found.

Does not invent allele, direction, PMID, population, or note values.

Example:
    uv run python scripts/curation/validate_longevity_allele_direction.py
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import polars as pl
import typer

from rogen_aging.config import find_repo_root
from rogen_aging.longevity.allele_harmonize import is_palindromic

REPO_ROOT = find_repo_root()
DEFAULT_CURATION = REPO_ROOT / "data" / "curation" / "longevity_allele_direction.csv"
DEFAULT_VARIANTS = REPO_ROOT / "data" / "processed" / "variants_47_input_v2.csv"
DEFAULT_PENDING = REPO_ROOT / "docs" / "PENDING_ALLELES.md"
DEFAULT_PROVENANCE = (
    REPO_ROOT / "data" / "curation" / "validate_longevity_allele_direction_provenance.json"
)
REQUIRED_SCHEMA = (
    "rsid",
    "gene",
    "longevity_associated_allele",
    "comparison",
    "source_pmid",
    "population",
    "notes",
)
# Completeness follows the GWAS-lookup docs: allele is required; PMID /
# population / notes are optional. ``comparison`` is optional metadata.
COMPLETE_REQUIRED = ("longevity_associated_allele",)
CONSTRAINED_COLUMNS = ("longevity_associated_allele", "comparison", "population")
PALINDROME_FLAG = "orientation cannot be inferred from frequency"
EXPECTED_N = 47


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_commit() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _repo_relative(path: Path) -> str:
    resolved = path.resolve()
    root = REPO_ROOT.resolve()
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def _cell(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() == "null":
        return ""
    return text


def load_curation(path: Path) -> pl.DataFrame:
    """Load the curation CSV and assert the ground-truth row count."""
    assert path.is_file(), f"MISSING INPUT: {path}"
    frame = pl.read_csv(path, schema_overrides={name: pl.Utf8 for name in REQUIRED_SCHEMA})
    missing = [name for name in REQUIRED_SCHEMA if name not in frame.columns]
    assert not missing, f"STOP: curation file missing columns: {missing}"
    assert (
        frame.height == EXPECTED_N
    ), f"ROW COUNT DISAGREES: computed={frame.height} ground_truth={EXPECTED_N}"
    return frame.select(list(REQUIRED_SCHEMA))


def load_panel(path: Path) -> pl.DataFrame:
    """Load panel alleles used only for palindrome flags."""
    assert path.is_file(), f"MISSING INPUT: {path}"
    frame = pl.read_csv(
        path,
        schema_overrides={
            "ref": pl.Utf8,
            "alt": pl.Utf8,
            "rsid": pl.Utf8,
            "gene_symbol": pl.Utf8,
        },
    )
    required = {"rsid", "gene_symbol", "ref", "alt"}
    missing = sorted(required - set(frame.columns))
    assert not missing, f"STOP: panel missing columns: {missing}"
    assert (
        frame.height == EXPECTED_N
    ), f"ROW COUNT DISAGREES: computed={frame.height} ground_truth={EXPECTED_N}"
    return frame


def completeness_report(frame: pl.DataFrame) -> tuple[list[str], list[str]]:
    """Return (complete_rsids, incomplete_rsids)."""
    complete: list[str] = []
    incomplete: list[str] = []
    for row in frame.iter_rows(named=True):
        rsid = _cell(row["rsid"])
        filled = all(_cell(row[name]) for name in COMPLETE_REQUIRED)
        if filled:
            complete.append(rsid)
        else:
            incomplete.append(rsid)
    return complete, incomplete


def distinct_values(frame: pl.DataFrame, columns: tuple[str, ...]) -> dict[str, list[str]]:
    """Collect sorted distinct non-empty values per column."""
    out: dict[str, list[str]] = {}
    for name in columns:
        values = sorted({_cell(value) for value in frame[name].to_list() if _cell(value)})
        out[name] = values
    return out


def documentation_paths(curation: Path, pending: Path) -> list[Path]:
    """Paths that may document allowed curation values."""
    candidates = [
        pending,
        curation.with_suffix(".md"),
        curation.parent / "README.md",
        curation.parent / "longevity_allele_direction.md",
        REPO_ROOT / "docs" / "LONGEVITY_GWAS_LOOKUP.md",
    ]
    return [path for path in candidates if path.is_file()]


def extract_allowed_values(doc_paths: list[Path]) -> dict[str, set[str]]:
    """Parse documented allow-lists. Returns {} when none are documented.

    Looks only for explicit blocks such as::

        allowed comparison values: a, b
        comparison: `x` | `y`

    Does not invent values from narrative text.
    """
    allowed: dict[str, set[str]] = {}
    patterns = [
        re.compile(
            r"allowed\s+(?P<column>longevity_associated_allele|comparison|population)"
            r"\s+values?\s*[:=]\s*(?P<body>.+)",
            re.IGNORECASE,
        ),
        re.compile(
            r"(?P<column>longevity_associated_allele|comparison|population)"
            r"\s+allowed\s+values?\s*[:=]\s*(?P<body>.+)",
            re.IGNORECASE,
        ),
    ]
    token_re = re.compile(r"[`'\"]([^`'\"]+)[`'\"]|([A-Za-z0-9_./+-]+)")

    for path in doc_paths:
        text = path.read_text(encoding="utf-8")
        for line in text.splitlines():
            for pattern in patterns:
                match = pattern.search(line)
                if match is None:
                    continue
                column = match.group("column").strip().lower()
                body = match.group("body").strip()
                tokens = {
                    token
                    for raw in token_re.findall(body)
                    for token in (raw[0] or raw[1],)
                    if token and token.lower() not in {"values", "value", "or", "and"}
                }
                if tokens:
                    allowed.setdefault(column, set()).update(tokens)
    return allowed


def palindrome_rows(panel: pl.DataFrame) -> list[dict[str, str]]:
    """Return panel rows whose ref/alt are A/T or C/G."""
    rows: list[dict[str, str]] = []
    for row in panel.iter_rows(named=True):
        ref = str(row["ref"] or "")
        alt = str(row["alt"] or "")
        if not is_palindromic(ref, alt):
            continue
        rows.append(
            {
                "rsid": str(row["rsid"]),
                "gene": str(row["gene_symbol"]),
                "ref": ref.upper(),
                "alt": alt.upper(),
                "flag": PALINDROME_FLAG,
            }
        )
    return rows


def write_provenance(
    path: Path,
    *,
    curation: Path,
    variants: Path,
    pending: Path,
    complete: list[str],
    incomplete: list[str],
    palindromes: list[dict[str, str]],
    distinct: dict[str, list[str]],
    doc_paths: list[Path],
    allowed: dict[str, set[str]],
    stop_reason: str | None,
) -> None:
    """Write tool versions, input hashes, git commit, and timestamp."""
    payload: dict[str, Any] = {
        "generated_at_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "git_commit": _git_commit(),
        "python": sys.version.split()[0],
        "polars": pl.__version__,
        "typer": typer.__version__,
        "expected_n_rows": EXPECTED_N,
        "complete_required_columns": list(COMPLETE_REQUIRED),
        "inputs": {
            "curation": {
                "path": _repo_relative(curation),
                "sha256": _sha256(curation),
            },
            "variants": {
                "path": _repo_relative(variants),
                "sha256": _sha256(variants),
            },
            "pending_alleles": {
                "path": _repo_relative(pending),
                "sha256": _sha256(pending),
            },
            "documentation_searched": [_repo_relative(path) for path in doc_paths],
        },
        "n_complete": len(complete),
        "n_incomplete": len(incomplete),
        "complete_rsids": complete,
        "incomplete_rsids": incomplete,
        "palindromic_variants": palindromes,
        "palindrome_flag": PALINDROME_FLAG,
        "distinct_values_found": distinct,
        "allowed_values_documented": {
            key: sorted(values) for key, values in sorted(allowed.items())
        },
        "stop_reason": stop_reason,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main(
    curation: Path = typer.Option(
        DEFAULT_CURATION,
        "--curation",
        help="longevity_allele_direction.csv to validate.",
    ),
    variants: Path = typer.Option(
        DEFAULT_VARIANTS,
        "--variants",
        help="Panel CSV used only to flag A/T and C/G palindromes.",
    ),
    pending: Path = typer.Option(
        DEFAULT_PENDING,
        "--pending",
        help="PENDING_ALLELES.md searched for documented allowed values.",
    ),
    provenance: Path = typer.Option(
        DEFAULT_PROVENANCE,
        "--provenance",
        help="Provenance JSON for this validation run.",
    ),
) -> None:
    """Report completeness, palindromes, and documented allow-lists."""
    assert pending.is_file(), f"MISSING INPUT: {pending}"
    frame = load_curation(curation)
    panel = load_panel(variants)
    complete, incomplete = completeness_report(frame)
    distinct = distinct_values(frame, CONSTRAINED_COLUMNS)
    palindromes = palindrome_rows(panel)
    doc_paths = documentation_paths(curation, pending)
    allowed = extract_allowed_values(doc_paths)

    typer.echo(f"complete={len(complete)}/{EXPECTED_N} incomplete={len(incomplete)}/{EXPECTED_N}")
    if incomplete:
        typer.echo("incomplete_rsids:")
        for rsid in incomplete:
            typer.echo(f"  {rsid}")
    else:
        typer.echo("incomplete_rsids: (none)")

    if palindromes:
        typer.echo(f"palindromic_variants={len(palindromes)} ({PALINDROME_FLAG})")
        for row in palindromes:
            typer.echo(
                f"  {row['rsid']} {row['gene']} {row['ref']}/{row['alt']}: {PALINDROME_FLAG}"
            )
    else:
        typer.echo("palindromic_variants=0")

    stop_reason: str | None = None
    if not allowed:
        stop_reason = (
            "No allowed values documented in PENDING_ALLELES.md or the "
            "curation file's own documentation."
        )
        write_provenance(
            provenance,
            curation=curation,
            variants=variants,
            pending=pending,
            complete=complete,
            incomplete=incomplete,
            palindromes=palindromes,
            distinct=distinct,
            doc_paths=doc_paths,
            allowed=allowed,
            stop_reason=stop_reason,
        )
        typer.echo(f"provenance={_repo_relative(provenance)}")
        typer.echo(stop_reason)
        typer.echo(f"documentation_searched={[_repo_relative(p) for p in doc_paths]}")
        typer.echo(f"distinct_values_found={json.dumps(distinct, sort_keys=True)}")
        assert False, f"{stop_reason} Distinct values found: {json.dumps(distinct, sort_keys=True)}"

    # Allow-list present: check non-empty constrained cells against it.
    violations: list[str] = []
    for row in frame.iter_rows(named=True):
        rsid = _cell(row["rsid"])
        for column in CONSTRAINED_COLUMNS:
            value = _cell(row[column])
            if not value:
                continue
            if column not in allowed:
                violations.append(f"{rsid} {column}={value!r} (no allow-list for this column)")
                continue
            if value not in allowed[column]:
                violations.append(f"{rsid} {column}={value!r} not in {sorted(allowed[column])}")

    write_provenance(
        provenance,
        curation=curation,
        variants=variants,
        pending=pending,
        complete=complete,
        incomplete=incomplete,
        palindromes=palindromes,
        distinct=distinct,
        doc_paths=doc_paths,
        allowed=allowed,
        stop_reason=None if not violations else "values outside documented allow-list",
    )
    typer.echo(f"provenance={_repo_relative(provenance)}")
    assert not violations, "STOP: values outside documented allow-list:\n" + "\n".join(violations)
    typer.echo("allow_list_check=passed")


if __name__ == "__main__":
    typer.run(main)
