"""Pre-flight validation for DNA methylation beta matrices used by aging clocks."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

from rogen_aging.clock.data import load_wide_table

SAMPLE_ID_CANDIDATES = ("sample_id", "Sample_ID", "SampleID", "IID", "id")


class MissingValuePolicy(StrEnum):
    """How to handle missing beta values at expected CpG sites."""

    REPORT = "report"
    FAIL = "fail"
    IMPUTE_COLUMN_MEAN = "impute_column_mean"
    DROP_SITES = "drop_sites"


@dataclass(frozen=True)
class ValidationIssue:
    """A single diagnostic finding from matrix preflight checks."""

    check: str
    severity: Literal["error", "warning", "info"]
    message: str
    count: int = 0
    examples: list[str] = field(default_factory=list)


@dataclass
class ValidationReport:
    """Summarized preflight diagnostics for a methylation matrix."""

    n_samples: int
    n_cpgs: int
    n_expected_cpgs: int
    n_missing_expected_cpgs: int
    n_extra_cpgs: int
    n_missing_beta_values: int
    n_out_of_range_beta_values: int
    n_samples_in_matrix_only: int
    n_samples_in_manifest_only: int
    n_duplicate_sample_ids: int
    missing_value_policy: str
    passed: bool
    issues: list[ValidationIssue] = field(default_factory=list)
    cleaned_n_cpgs: int | None = None

    def to_dict(self) -> dict[str, object]:
        """Serialize the report to a JSON-friendly dict."""
        payload = asdict(self)
        return payload

    def format_log(self) -> str:
        """Return a clean, human-readable diagnostic summary."""
        status = "PASSED" if self.passed else "FAILED"
        lines = [
            "=== Methylation matrix preflight ===",
            f"Status: {status}",
            f"Samples: {self.n_samples}",
            f"CpG columns in matrix: {self.n_cpgs}",
            f"Expected CpG sites: {self.n_expected_cpgs}",
            f"Missing expected CpGs: {self.n_missing_expected_cpgs}",
            f"Extra (unexpected) CpGs: {self.n_extra_cpgs}",
            f"Missing beta values (cells): {self.n_missing_beta_values}",
            f"Out-of-range beta values (not in (0, 1)): {self.n_out_of_range_beta_values}",
            f"Sample IDs only in matrix: {self.n_samples_in_matrix_only}",
            f"Sample IDs only in manifest: {self.n_samples_in_manifest_only}",
            f"Duplicate sample IDs: {self.n_duplicate_sample_ids}",
            f"Missing-value policy: {self.missing_value_policy}",
        ]
        if self.cleaned_n_cpgs is not None:
            lines.append(f"CpGs after cleaning: {self.cleaned_n_cpgs}")
        if self.issues:
            lines.append("")
            lines.append("--- Findings ---")
            for issue in self.issues:
                prefix = issue.severity.upper()
                count_part = f" (n={issue.count})" if issue.count else ""
                lines.append(f"[{prefix}] {issue.check}: {issue.message}{count_part}")
                if issue.examples:
                    shown = ", ".join(issue.examples[:8])
                    more = "" if len(issue.examples) <= 8 else f" … (+{len(issue.examples) - 8})"
                    lines.append(f"         examples: {shown}{more}")
        return "\n".join(lines)


def _resolve_sample_id_column(df: pd.DataFrame, sample_id_col: str | None) -> str:
    if sample_id_col is not None:
        if sample_id_col not in df.columns:
            raise ValueError(f"Sample ID column '{sample_id_col}' not found in table.")
        return sample_id_col
    for candidate in SAMPLE_ID_CANDIDATES:
        if candidate in df.columns:
            return candidate
    raise ValueError(
        "No sample ID column found. Expected one of: "
        + ", ".join(SAMPLE_ID_CANDIDATES)
        + ". Pass sample_id_col explicitly if needed."
    )


def _cpg_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if str(c).startswith("cg")]


def _load_table(path: Path) -> pd.DataFrame:
    return load_wide_table(path)


def _load_expected_cpgs(path: Path | None, matrix_cpgs: list[str]) -> list[str]:
    if path is None:
        return list(matrix_cpgs)
    if not path.is_file():
        raise FileNotFoundError(path)
    text = path.read_text(encoding="utf-8")
    names: list[str] = []
    for line in text.splitlines():
        token = line.strip().split(",")[0].strip().strip('"')
        if not token or token.lower() in {"cpg", "probe", "ilmnid", "id_ref", "feature"}:
            continue
        names.append(token)
    if not names:
        raise ValueError(f"No CpG identifiers found in expected-CpG file: {path}")
    # Preserve order, drop duplicates
    seen: set[str] = set()
    ordered: list[str] = []
    for name in names:
        if name not in seen:
            seen.add(name)
            ordered.append(name)
    return ordered


def _sample_ids(series: pd.Series) -> list[str]:
    return series.astype(str).str.strip().tolist()


def _apply_missing_policy(
    beta: pd.DataFrame,
    policy: MissingValuePolicy,
    issues: list[ValidationIssue],
    *,
    examples: list[str] | None = None,
) -> pd.DataFrame:
    missing_cells = int(beta.isna().sum().sum())
    if missing_cells == 0:
        return beta

    example_list = examples or []

    if policy is MissingValuePolicy.FAIL:
        issues.append(
            ValidationIssue(
                check="missing_values",
                severity="error",
                message="Missing beta values present and policy is 'fail'.",
                count=missing_cells,
                examples=example_list,
            )
        )
        return beta

    if policy is MissingValuePolicy.REPORT:
        issues.append(
            ValidationIssue(
                check="missing_values",
                severity="warning",
                message="Missing beta values reported; matrix left unchanged.",
                count=missing_cells,
                examples=example_list,
            )
        )
        return beta

    if policy is MissingValuePolicy.DROP_SITES:
        keep = beta.columns[beta.notna().all()]
        dropped = [c for c in beta.columns if c not in set(keep)]
        issues.append(
            ValidationIssue(
                check="missing_values",
                severity="warning",
                message="Dropped CpG sites with any missing values.",
                count=len(dropped),
                examples=[str(c) for c in dropped[:20]],
            )
        )
        return beta.loc[:, keep].copy()

    if policy is MissingValuePolicy.IMPUTE_COLUMN_MEAN:
        filled = beta.copy()
        col_means = filled.mean(axis=0, skipna=True)
        global_mean = float(col_means.mean()) if col_means.notna().any() else 0.5
        if not pd.notna(global_mean):
            global_mean = 0.5
        for col in filled.columns:
            fill = col_means[col]
            if pd.isna(fill):
                fill = global_mean
            filled[col] = filled[col].fillna(float(fill))
        issues.append(
            ValidationIssue(
                check="missing_values",
                severity="warning",
                message="Imputed missing beta values with per-CpG column means.",
                count=missing_cells,
                examples=example_list,
            )
        )
        return filled

    raise ValueError(f"Unsupported missing-value policy: {policy}")


def validate_methylation_matrix(
    matrix: pd.DataFrame | Path,
    metadata: pd.DataFrame | Path,
    *,
    expected_cpgs: Path | list[str] | None = None,
    sample_id_col: str | None = None,
    metadata_sample_id_col: str | None = None,
    missing_policy: MissingValuePolicy | str = MissingValuePolicy.REPORT,
    require_strict_beta: bool = True,
) -> tuple[ValidationReport, pd.DataFrame]:
    """Validate a methylation beta matrix against a metadata manifest.

    Checks:
      1. Missing values at expected CpG sites (report / fail / impute / drop).
      2. Beta values strictly inside ``(0, 1)`` when ``require_strict_beta``.
      3. Sample IDs match the metadata manifest (set equality after de-dup check).

    Args:
        matrix: Wide beta table or path (``.parquet`` / ``.csv`` / ``.tsv``).
        metadata: Manifest table or path with sample IDs.
        expected_cpgs: Optional CpG list (path or names). Defaults to matrix ``cg*``.
        sample_id_col: Sample ID column in the matrix (auto-detected if omitted).
        metadata_sample_id_col: Sample ID column in metadata (auto-detected).
        missing_policy: How to handle missing betas at expected sites.
        require_strict_beta: If true, values must satisfy ``0 < beta < 1``.

    Returns:
        A pair ``(report, cleaned_matrix)`` where ``cleaned_matrix`` retains
        ``sample_id`` plus CpG columns after applying the missing-value policy.

    Raises:
        ValueError: If required columns are missing or tables are empty.
        FileNotFoundError: If a provided path does not exist.
    """
    policy = MissingValuePolicy(missing_policy)
    issues: list[ValidationIssue] = []

    matrix_df = _load_table(matrix) if isinstance(matrix, Path) else matrix.copy()
    meta_df = _load_table(metadata) if isinstance(metadata, Path) else metadata.copy()

    if matrix_df.empty:
        raise ValueError("Methylation matrix is empty.")
    if meta_df.empty:
        raise ValueError("Metadata manifest is empty.")

    matrix_id_col = _resolve_sample_id_column(matrix_df, sample_id_col)
    meta_id_col = _resolve_sample_id_column(meta_df, metadata_sample_id_col)

    matrix_ids = _sample_ids(matrix_df[matrix_id_col])
    meta_ids = _sample_ids(meta_df[meta_id_col])

    matrix_id_counts = pd.Series(matrix_ids).value_counts()
    meta_id_counts = pd.Series(meta_ids).value_counts()
    dup_matrix = sorted(matrix_id_counts[matrix_id_counts > 1].index.tolist())
    dup_meta = sorted(meta_id_counts[meta_id_counts > 1].index.tolist())
    n_duplicate_sample_ids = len(set(dup_matrix) | set(dup_meta))
    if dup_matrix:
        issues.append(
            ValidationIssue(
                check="sample_ids",
                severity="error",
                message="Duplicate sample IDs in methylation matrix.",
                count=len(dup_matrix),
                examples=dup_matrix[:20],
            )
        )
    if dup_meta:
        issues.append(
            ValidationIssue(
                check="sample_ids",
                severity="error",
                message="Duplicate sample IDs in metadata manifest.",
                count=len(dup_meta),
                examples=dup_meta[:20],
            )
        )

    matrix_id_set = set(matrix_ids)
    meta_id_set = set(meta_ids)
    only_matrix = sorted(matrix_id_set - meta_id_set)
    only_meta = sorted(meta_id_set - matrix_id_set)
    if only_matrix:
        issues.append(
            ValidationIssue(
                check="sample_ids",
                severity="error",
                message="Sample IDs present in matrix but missing from metadata manifest.",
                count=len(only_matrix),
                examples=only_matrix[:20],
            )
        )
    if only_meta:
        issues.append(
            ValidationIssue(
                check="sample_ids",
                severity="error",
                message="Sample IDs present in metadata manifest but missing from matrix.",
                count=len(only_meta),
                examples=only_meta[:20],
            )
        )

    matrix_cpgs = _cpg_columns(matrix_df)
    if not matrix_cpgs:
        raise ValueError("No CpG columns starting with 'cg' found in the methylation matrix.")

    if isinstance(expected_cpgs, list):
        expected = list(dict.fromkeys(expected_cpgs))
    else:
        expected = _load_expected_cpgs(expected_cpgs, matrix_cpgs)

    expected_set = set(expected)
    matrix_cpg_set = set(matrix_cpgs)
    missing_expected = [c for c in expected if c not in matrix_cpg_set]
    extra_cpgs = [c for c in matrix_cpgs if c not in expected_set]

    if missing_expected:
        issues.append(
            ValidationIssue(
                check="expected_cpgs",
                severity="error",
                message="Expected CpG sites absent from the methylation matrix.",
                count=len(missing_expected),
                examples=missing_expected[:20],
            )
        )
    if extra_cpgs and expected_cpgs is not None:
        issues.append(
            ValidationIssue(
                check="expected_cpgs",
                severity="info",
                message="Matrix contains CpG sites not listed in the expected-CpG set.",
                count=len(extra_cpgs),
                examples=extra_cpgs[:20],
            )
        )

    present_expected = [c for c in expected if c in matrix_cpg_set]
    beta = matrix_df.reindex(columns=present_expected).apply(pd.to_numeric, errors="coerce")
    n_missing_beta = int(beta.isna().sum().sum())
    missing_examples: list[str] = []
    if n_missing_beta:
        site_missing = beta.isna().sum()
        worst = site_missing[site_missing > 0].sort_values(ascending=False)
        missing_examples = [f"{idx}:{int(val)}" for idx, val in worst.head(12).items()]

    beta = _apply_missing_policy(beta, policy, issues, examples=missing_examples)

    finite = beta.to_numpy(dtype=float, copy=False)
    # NaN is missing, not out-of-range.
    mask_finite = pd.notna(beta).to_numpy()
    if require_strict_beta:
        out_of_range_mask = mask_finite & ~((finite > 0.0) & (finite < 1.0))
    else:
        out_of_range_mask = mask_finite & ~((finite >= 0.0) & (finite <= 1.0))
    n_out_of_range = int(out_of_range_mask.sum())
    if n_out_of_range:
        rr, cc = np.nonzero(out_of_range_mask)
        sample_labels = matrix_df[matrix_id_col].astype(str).to_numpy()
        examples = [
            f"{sample_labels[int(r_i)]}:{beta.columns[int(c_i)]}={finite[int(r_i), int(c_i)]:.6g}"
            for r_i, c_i in zip(rr[:12], cc[:12], strict=False)
        ]
        issues.append(
            ValidationIssue(
                check="beta_range",
                severity="error",
                message=(
                    "Beta values must fall strictly between 0 and 1 (exclusive)."
                    if require_strict_beta
                    else "Beta values must fall within [0, 1]."
                ),
                count=n_out_of_range,
                examples=examples,
            )
        )

    cleaned = pd.DataFrame({matrix_id_col: matrix_df[matrix_id_col].astype(str).str.strip()})
    cleaned = pd.concat([cleaned, beta.reset_index(drop=True)], axis=1)

    passed = not any(issue.severity == "error" for issue in issues)
    report = ValidationReport(
        n_samples=len(matrix_df),
        n_cpgs=len(matrix_cpgs),
        n_expected_cpgs=len(expected),
        n_missing_expected_cpgs=len(missing_expected),
        n_extra_cpgs=len(extra_cpgs) if expected_cpgs is not None else 0,
        n_missing_beta_values=n_missing_beta,
        n_out_of_range_beta_values=n_out_of_range,
        n_samples_in_matrix_only=len(only_matrix),
        n_samples_in_manifest_only=len(only_meta),
        n_duplicate_sample_ids=n_duplicate_sample_ids,
        missing_value_policy=policy.value,
        passed=passed,
        issues=issues,
        cleaned_n_cpgs=len(beta.columns),
    )
    return report, cleaned


def write_validation_outputs(
    report: ValidationReport,
    cleaned: pd.DataFrame,
    *,
    log_path: Path | None = None,
    report_json_path: Path | None = None,
    cleaned_matrix_path: Path | None = None,
) -> None:
    """Write diagnostic log, optional JSON report, and optional cleaned matrix."""
    log_text = report.format_log()
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(log_text + "\n", encoding="utf-8")
    if report_json_path is not None:
        report_json_path.parent.mkdir(parents=True, exist_ok=True)
        report_json_path.write_text(
            json.dumps(report.to_dict(), indent=2) + "\n",
            encoding="utf-8",
        )
    if cleaned_matrix_path is not None:
        cleaned_matrix_path.parent.mkdir(parents=True, exist_ok=True)
        suffix = cleaned_matrix_path.suffix.lower()
        if suffix == ".parquet":
            cleaned.to_parquet(cleaned_matrix_path, index=False)
        elif suffix == ".tsv":
            cleaned.to_csv(cleaned_matrix_path, sep="\t", index=False)
        else:
            cleaned.to_csv(cleaned_matrix_path, index=False)
