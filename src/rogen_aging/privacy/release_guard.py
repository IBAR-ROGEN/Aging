"""Release checks for aggregate tables leaving the processing zone.

Checks 1 and 2 (individual-level data and direct identifiers) are never
rewritten. A failure blocks release. Checks 3 and 4 suppress small cells
on the sanitized copy only.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import polars as pl
from omegaconf import DictConfig, OmegaConf

CheckStatus = Literal["PASS", "FAIL"]
TableKind = Literal["allele_counts", "strata", "summary"]

_TABLE_KINDS: frozenset[str] = frozenset({"allele_counts", "strata", "summary"})
_AC_NAMES: frozenset[str] = frozenset({"ac", "allele_count"})
_AN_NAMES: frozenset[str] = frozenset({"an", "allele_number"})
_AF_NAMES: frozenset[str] = frozenset({"af", "allele_frequency"})
_TOTAL_NAMES: frozenset[str] = frozenset({"total", "row_total", "n_total"})
_NAME_COLUMNS: frozenset[str] = frozenset(
    {
        "name",
        "full_name",
        "first_name",
        "last_name",
        "surname",
        "forename",
        "given_name",
        "patient_name",
        "family_name",
    }
)
_POSTCODE_COLUMNS: frozenset[str] = frozenset(
    {"postcode", "post_code", "postal_code", "zip", "zip_code", "zipcode"}
)
_DOB_COLUMNS: frozenset[str] = frozenset(
    {"date_of_birth", "birth_date", "birthdate", "dob", "d_o_b"}
)

_INDIVIDUAL_COLUMN = re.compile(r"(?:^|_)(?:sample_?id|eid|participant|pseudonym|iid|fid)(?:_|$)")
_PSEUDONYM = re.compile(r"(?<![0-9A-Fa-f])RO-[0-9A-Fa-f]{12}(?![0-9A-Fa-f])")
_FULL_DATE_VALUE = re.compile(
    r"^(?:\d{4}-\d{2}-\d{2}(?:[T ]\d{2}:\d{2}:\d{2}(?:\.\d+)?Z?)?|\d{2}/\d{2}/\d{4})$"
)
_POSTCODE_VALUE = re.compile(r"^[A-Z]{1,2}\d[A-Z\d]?\s*\d[A-Z]{2}$", re.IGNORECASE)


class ReleaseGuardError(Exception):
    """The release guard could not run on this table."""


@dataclass(frozen=True)
class CheckResult:
    """One PASS or FAIL line in the release report."""

    check: str
    status: CheckStatus
    reason: str

    def as_dict(self) -> dict[str, str]:
        return {"check": self.check, "status": self.status, "reason": self.reason}


@dataclass(frozen=True)
class ReleaseGuardResult:
    """Checks plus the sanitized frame when release is not blocked."""

    table_kind: TableKind
    min_cell_count: int
    n_samples: int | None
    checks: tuple[CheckResult, ...]
    sanitized: pl.DataFrame | None

    @property
    def release_blocked(self) -> bool:
        """True when any check failed.

        Checks 1 and 2 always block, and they are not repaired. Checks 3 and 4
        fail only when suppression cannot be applied; that also blocks release.
        """
        return any(item.status == "FAIL" for item in self.checks)

    def report(self, *, table: str, sanitized_written: bool) -> dict[str, object]:
        return {
            "table": table,
            "table_kind": self.table_kind,
            "min_cell_count": self.min_cell_count,
            "n_samples": self.n_samples,
            "release_blocked": self.release_blocked,
            "sanitized_written": sanitized_written,
            "checks": [item.as_dict() for item in self.checks],
        }


def _normalize(name: str) -> str:
    return re.sub(r"[\s\-]+", "_", name.strip().lower())


def _cell_strings(frame: pl.DataFrame) -> dict[str, list[str]]:
    found: dict[str, list[str]] = {}
    for name in frame.columns:
        raw = frame.get_column(name).cast(pl.Utf8, strict=False).to_list()
        texts: list[str] = []
        for value in raw:
            if value is None:
                continue
            text = str(value).strip()
            if text:
                texts.append(text)
        found[name] = texts
    return found


def _match_column(columns: list[str], names: frozenset[str]) -> str | None:
    for column in columns:
        if _normalize(column) in names:
            return column
    return None


def check_individual_level(frame: pl.DataFrame, *, n_samples: int | None) -> CheckResult:
    """Fail on individual identifiers, RO- pseudonyms, or one row per sample."""
    reasons: list[str] = []
    named = [
        column
        for column in frame.columns
        if _INDIVIDUAL_COLUMN.search(_normalize(column)) is not None
    ]
    if named:
        joined = ", ".join(named)
        reasons.append(f"column name matches an individual-level field: {joined}")

    texts = _cell_strings(frame)
    hidden = [
        column
        for column, values in texts.items()
        if any(_PSEUDONYM.search(value) for value in values)
    ]
    if hidden:
        joined = ", ".join(hidden)
        reasons.append(f"cell matches pseudonym pattern RO- plus 12 hex characters in: {joined}")

    if n_samples is not None and frame.height == n_samples:
        reasons.append(
            f"row count {frame.height} equals --n-samples {n_samples}; table has one row per sample"
        )

    if reasons:
        return CheckResult("individual_level", "FAIL", "; ".join(reasons))
    if n_samples is None:
        reason = "no individual-level column and no RO- pseudonym; --n-samples was not supplied"
    else:
        reason = (
            "no individual-level column, no RO- pseudonym, "
            f"and row count {frame.height} does not equal --n-samples {n_samples}"
        )
    return CheckResult("individual_level", "PASS", reason)


def _is_dob_column(name: str) -> bool:
    return _normalize(name) in _DOB_COLUMNS


def _is_date_column(name: str, values: list[str]) -> bool:
    normalized = _normalize(name)
    if normalized == "date" or normalized.endswith("_date") or normalized == "full_date":
        return True
    if not values:
        return False
    return all(_FULL_DATE_VALUE.fullmatch(value) is not None for value in values)


def _is_postcode_column(name: str, values: list[str]) -> bool:
    if _normalize(name) in _POSTCODE_COLUMNS:
        return True
    if not values:
        return False
    return all(_POSTCODE_VALUE.fullmatch(value) is not None for value in values)


def _is_name_column(name: str) -> bool:
    return _normalize(name) in _NAME_COLUMNS


def check_direct_identifiers(frame: pl.DataFrame) -> CheckResult:
    """Fail on columns that look like dates of birth, full dates, postcodes, or names."""
    texts = _cell_strings(frame)
    reasons: list[str] = []
    for column in frame.columns:
        values = texts[column]
        if _is_dob_column(column):
            reasons.append(f"column {column} looks like a date of birth")
        elif _is_date_column(column, values):
            reasons.append(f"column {column} looks like full dates")
        elif _is_postcode_column(column, values):
            reasons.append(f"column {column} looks like postcodes")
        elif _is_name_column(column):
            reasons.append(f"column {column} looks like names")
    if reasons:
        return CheckResult("direct_identifiers", "FAIL", "; ".join(reasons))
    return CheckResult(
        "direct_identifiers",
        "PASS",
        "no date-of-birth, full-date, postcode, or name column",
    )


def _count_columns(frame: pl.DataFrame) -> tuple[list[str], list[str]]:
    counts: list[str] = []
    totals: list[str] = []
    for name, dtype in zip(frame.columns, frame.dtypes, strict=True):
        if not dtype.is_numeric():
            continue
        if _normalize(name) in _TOTAL_NAMES:
            totals.append(name)
        else:
            counts.append(name)
    return counts, totals


def suppress_allele_counts(
    frame: pl.DataFrame,
    *,
    min_cell_count: int,
) -> tuple[pl.DataFrame, CheckResult]:
    """Drop rows with a small allele count or a small complementary count."""
    ac_name = _match_column(frame.columns, _AC_NAMES)
    an_name = _match_column(frame.columns, _AN_NAMES)
    if ac_name is None or an_name is None:
        missing = []
        if ac_name is None:
            missing.append("AC or allele_count")
        if an_name is None:
            missing.append("AN or allele_number")
        reason = "allele_counts table is missing " + " and ".join(missing)
        return frame, CheckResult("allele_count_cells", "FAIL", reason)

    small = (pl.col(ac_name) < min_cell_count) | (
        (pl.col(an_name) - pl.col(ac_name)) < min_cell_count
    )
    small = small | pl.col(ac_name).is_null() | pl.col(an_name).is_null()
    kept = frame.filter(~small)
    n_suppressed = frame.height - kept.height
    af_names = [column for column in kept.columns if _normalize(column) in _AF_NAMES]
    if af_names:
        kept = kept.with_columns([pl.col(name).round(3) for name in af_names])
        rounded = "allele frequencies rounded to 3 decimals"
    else:
        rounded = "no allele-frequency column to round"
    reason = (
        f"suppressed {n_suppressed} row(s) where AC < {min_cell_count} "
        f"or AN - AC < {min_cell_count}; {rounded}"
    )
    return kept, CheckResult("allele_count_cells", "PASS", reason)


def suppress_strata(
    frame: pl.DataFrame,
    *,
    min_cell_count: int,
) -> tuple[pl.DataFrame, CheckResult]:
    """Suppress counts below k, then one more cell when the row total would recover it."""
    count_cols, total_cols = _count_columns(frame)
    if not count_cols:
        return frame, CheckResult(
            "strata_cells",
            "FAIL",
            "strata table has no numeric phenotype-group count columns",
        )

    records = frame.select([*count_cols, *total_cols]).to_dicts()
    n_primary = 0
    n_extra = 0
    for record in records:
        small = [
            column
            for column in count_cols
            if record[column] is not None and float(record[column]) < min_cell_count
        ]
        n_primary += len(small)
        extra: str | None = None
        if len(small) == 1:
            remaining = [
                column
                for column in count_cols
                if column not in small and record[column] is not None
            ]
            if remaining:
                extra = min(
                    remaining,
                    key=lambda column: (float(record[column]), count_cols.index(column)),
                )
            elif total_cols:
                extra = total_cols[0]
        if extra is not None:
            n_extra += 1
            small.append(extra)
        for column in small:
            record[column] = None

    updated = frame.with_columns(
        [
            pl.Series(column, [record[column] for record in records], dtype=frame.schema[column])
            for column in (*count_cols, *total_cols)
        ]
    )
    reason = (
        f"suppressed {n_primary} cell(s) below {min_cell_count}; "
        f"complementary suppression hid {n_extra} additional cell(s) "
        "whose row total would have recovered a single hidden count"
    )
    return updated, CheckResult("strata_cells", "PASS", reason)


def run_release_guard(
    frame: pl.DataFrame,
    *,
    table_kind: str,
    min_cell_count: int,
    n_samples: int | None = None,
) -> ReleaseGuardResult:
    """Run every release check. Sanitized output is omitted when release is blocked."""
    if table_kind not in _TABLE_KINDS:
        allowed = ", ".join(sorted(_TABLE_KINDS))
        raise ReleaseGuardError(f"table kind must be one of: {allowed}")
    if min_cell_count < 1:
        raise ReleaseGuardError("min_cell_count must be at least 1")
    if table_kind == "allele_counts":
        kind: TableKind = "allele_counts"
    elif table_kind == "strata":
        kind = "strata"
    else:
        kind = "summary"

    individual = check_individual_level(frame, n_samples=n_samples)
    identifiers = check_direct_identifiers(frame)
    sanitized: pl.DataFrame | None
    if kind == "allele_counts":
        sanitized, allele = suppress_allele_counts(frame, min_cell_count=min_cell_count)
        strata = CheckResult(
            "strata_cells",
            "PASS",
            "not applied: table kind is allele_counts",
        )
    elif kind == "strata":
        sanitized, strata = suppress_strata(frame, min_cell_count=min_cell_count)
        allele = CheckResult(
            "allele_count_cells",
            "PASS",
            "not applied: table kind is strata",
        )
    else:
        sanitized = frame
        allele = CheckResult(
            "allele_count_cells",
            "PASS",
            "not applied: table kind is summary",
        )
        strata = CheckResult(
            "strata_cells",
            "PASS",
            "not applied: table kind is summary",
        )

    checks = (individual, identifiers, allele, strata)
    if any(item.status == "FAIL" for item in checks):
        sanitized = None
    return ReleaseGuardResult(
        table_kind=kind,
        min_cell_count=min_cell_count,
        n_samples=n_samples,
        checks=checks,
        sanitized=sanitized,
    )


def load_table(path: Path) -> pl.DataFrame:
    """Read a CSV, TSV, or Parquet table."""
    if not path.is_file():
        raise ReleaseGuardError(f"table not found: {path}")
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pl.read_parquet(path)
    if suffix == ".tsv":
        return pl.read_csv(path, separator="\t")
    if suffix == ".csv":
        return pl.read_csv(path)
    raise ReleaseGuardError(f"unsupported table type: {path.suffix}")


def write_table(frame: pl.DataFrame, path: Path) -> None:
    """Write a sanitized CSV, TSV, or Parquet table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        frame.write_parquet(path)
        return
    if suffix == ".tsv":
        frame.write_csv(path, separator="\t")
        return
    if suffix == ".csv":
        frame.write_csv(path)
        return
    raise ReleaseGuardError(f"unsupported sanitized table type: {path.suffix}")


def write_report(payload: dict[str, object], path: Path) -> None:
    """Write the JSON release report."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def min_cell_count_from_config(cfg: object) -> int:
    """Read ``release_guard.min_cell_count`` from a loaded config."""
    if not isinstance(cfg, DictConfig):
        raise ReleaseGuardError("config must be an OmegaConf DictConfig")
    value = OmegaConf.select(cfg, "release_guard.min_cell_count")
    if value is None:
        raise ReleaseGuardError("release_guard.min_cell_count is missing from config")
    count = int(value)
    if count < 1:
        raise ReleaseGuardError("release_guard.min_cell_count must be at least 1")
    return count
