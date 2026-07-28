"""Synthetic UKB phenotype–genotype join and LA-SNP association scan.

Activity **2.1.11.1** — architecture validation on **synthetic** mock RAP output
only (from ``scripts/ukb/mock_rap_folder.py``). No real UK Biobank data; outputs
are for pipeline QA, **not** biological conclusions.

This package (``rogen_aging.ukb_integration``) is distinct from
``rogen_aging.integrative``, which performs offline variant×tissue×phenotype
risk scoring on curated annotation tables.
"""

from __future__ import annotations

import logging
import math
from datetime import UTC, datetime
from pathlib import Path
from typing import Final

import cyvcf2
import numpy as np
import polars as pl
from scipy.stats import fisher_exact, norm

LOG = logging.getLogger(__name__)

ACTIVITY_ID: Final[str] = "2.1.11.1"
SYNTHETIC_DISCLAIMER: Final[str] = (
    "Synthetic-data validation only (Activity 2.1.11.1); do not interpret biologically."
)
DEFAULT_AUDIT_LOG: Final[Path] = Path("outputs/logs/ukb_integration_audit.log")
MAX_JOIN_DROP_RATE: Final[float] = 0.01
_EID_AUDIT_SAMPLE_LIMIT: Final[int] = 50


class JoinDropRateError(ValueError):
    """Raised when phenotype–genotype ``eid`` join drop rate exceeds the allowed threshold."""

LA_SNP_ASSOC_COLUMNS: Final[tuple[str, ...]] = (
    "rsID",
    "OR",
    "CI_low",
    "CI_high",
    "p_value",
    "n",
)

PARENTAL_LONGEVITY_OUT = "assoc_la_snp_parental_longevity.csv"
AD_OUT = "assoc_la_snp_ad.csv"

PhenotypeColumn = str


def load_phenotype_table(path: Path) -> pl.DataFrame:
    """Load mock UKB phenotype CSV (``eid`` + v2 fields) with safety comment lines skipped.

    Args:
        path: Path to the phenotype CSV (``#``-prefixed lines are ignored).

    Returns:
        Polars frame with at least an ``eid`` column and mock v2 phenotype fields.

    Raises:
        FileNotFoundError: If ``path`` does not exist or is not a file.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Phenotype table not found: {path.resolve()}")
    return pl.read_csv(path, comment_prefix="#")


def _alt_dosage_from_gt_type(gt_type: int) -> int | None:
    """Map cyvcf2 ``gt_types`` to alt-allele dosage 0/1/2.

    cyvcf2 encoding: ``0=HOM_REF``, ``1=HET``, ``2=UNKNOWN``, ``3=HOM_ALT``.

    Args:
        gt_type: Integer genotype class from ``variant.gt_types``.

    Returns:
        Alt-allele dosage ``0``, ``1``, or ``2``, or ``None`` for unknown/missing.
    """
    if gt_type == 0:
        return 0
    if gt_type == 1:
        return 1
    if gt_type == 3:
        return 2
    return None


def load_genotype_matrix_from_vcf(path: Path) -> pl.DataFrame:
    """Load LA-SNP VCF into a wide matrix: one row per ``eid``, one column per rsID (0/1/2).

    Args:
        path: Path to a VCF whose sample IDs match phenotype ``eid`` values.

    Returns:
        Wide Polars frame with ``eid`` and one integer dosage column per variant
        (missing/unknown genotypes as null).

    Raises:
        FileNotFoundError: If ``path`` does not exist or is not a file.
        ValueError: If the VCF has no samples or no variant records.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Genotype VCF not found: {path.resolve()}")

    vcf = cyvcf2.VCF(str(path))
    sample_ids = list(vcf.samples)
    if not sample_ids:
        raise ValueError(f"VCF has no samples: {path}")

    dosage_by_snp: dict[str, list[int | None]] = {}
    for variant in vcf:
        rs_id = variant.ID
        if not rs_id:
            chrom = variant.CHROM
            pos = variant.POS
            rs_id = f"{chrom}:{pos}"
        dosages: list[int | None] = [
            _alt_dosage_from_gt_type(int(gt_type)) for gt_type in variant.gt_types
        ]
        dosage_by_snp[str(rs_id)] = dosages

    if not dosage_by_snp:
        raise ValueError(f"VCF contains no variant records: {path}")

    return pl.DataFrame({"eid": sample_ids, **dosage_by_snp})


def _eid_dtype_name(frame: pl.DataFrame) -> str:
    """Return the Polars dtype name of the ``eid`` column."""
    return str(frame.schema["eid"])


def _normalize_eid_frame(frame: pl.DataFrame) -> pl.DataFrame:
    """Cast ``eid`` to Utf8 and strip surrounding whitespace for stable joins."""
    return frame.with_columns(
        pl.col("eid").cast(pl.Utf8).str.strip_chars().alias("eid")
    )


def _sorted_unique_eids(frame: pl.DataFrame) -> list[str]:
    """Return sorted unique non-null ``eid`` values as strings."""
    return sorted(frame.get_column("eid").drop_nulls().unique().to_list())


def _format_eid_sample(eids: list[str], *, limit: int = _EID_AUDIT_SAMPLE_LIMIT) -> str:
    """Format a capped sample of EIDs for the audit log."""
    if not eids:
        return "(none)"
    shown = eids[:limit]
    suffix = f" … (+{len(eids) - limit} more)" if len(eids) > limit else ""
    return ", ".join(shown) + suffix


def write_join_audit_log(
    path: Path,
    *,
    pheno_n: int,
    geno_n: int,
    matched_n: int,
    pheno_only: list[str],
    geno_only: list[str],
    pheno_dtype: str,
    geno_dtype: str,
    schema_mismatch: bool,
    schema_normalized: bool,
    drop_rate: float,
    max_drop_rate: float,
    halted: bool,
) -> None:
    """Write phenotype–genotype ``eid`` join diagnostics to ``path``.

    Args:
        path: Destination audit log path; parent directories are created as needed.
        pheno_n: Unique phenotype ``eid`` count before join.
        geno_n: Unique genotype sample-ID count before join.
        matched_n: Size of the ``eid`` intersection after normalization.
        pheno_only: Phenotype ``eid`` values absent from genotypes.
        geno_only: Genotype sample IDs absent from phenotypes.
        pheno_dtype: Original phenotype ``eid`` Polars dtype name.
        geno_dtype: Original genotype ``eid`` Polars dtype name.
        schema_mismatch: Whether original ``eid`` dtypes differed.
        schema_normalized: Whether both sides were cast to Utf8 before joining.
        drop_rate: Fraction of union IDs that failed to match.
        max_drop_rate: Threshold above which the pipeline must halt.
        halted: Whether the join was aborted due to excess drop rate.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = [
        f"[{stamp}] Activity {ACTIVITY_ID} UKB phenotype–genotype eid join audit",
        f"phenotype_unique_eids={pheno_n}",
        f"genotype_unique_eids={geno_n}",
        f"matched_eids={matched_n}",
        f"pheno_only_dropped={len(pheno_only)}",
        f"geno_only_dropped={len(geno_only)}",
        f"drop_rate={drop_rate:.6f}",
        f"max_drop_rate={max_drop_rate:.6f}",
        f"eid_schema_phenotype={pheno_dtype}",
        f"eid_schema_genotype={geno_dtype}",
        f"eid_schema_mismatch={schema_mismatch}",
        f"eid_schema_normalized_to_utf8={schema_normalized}",
        f"halted={halted}",
        f"pheno_only_eids_sample={_format_eid_sample(pheno_only)}",
        f"geno_only_eids_sample={_format_eid_sample(geno_only)}",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    LOG.info("Wrote eid join audit log → %s (drop_rate=%.4f)", path, drop_rate)


def join_phenotypes_genotypes(
    phenotypes: pl.DataFrame,
    genotypes: pl.DataFrame,
    *,
    audit_log: Path | None = DEFAULT_AUDIT_LOG,
    max_drop_rate: float = MAX_JOIN_DROP_RATE,
) -> pl.DataFrame:
    """Inner-join phenotype table and genotype matrix on ``eid`` (one row per participant).

    Normalizes ``eid`` to stripped Utf8 on both sides so Int/Utf8 schema mismatches
    do not silently zero the join. When ``audit_log`` is set, dropped IDs and schema
    diagnostics are written there. The pipeline halts if the union drop rate exceeds
    ``max_drop_rate`` (default 1%).

    Args:
        phenotypes: Phenotype frame with an ``eid`` column.
        genotypes: Genotype matrix with an ``eid`` column and SNP dosage columns.
        audit_log: Path for the join audit log, or ``None`` to skip file output.
        max_drop_rate: Maximum allowed fraction of unmatched union ``eid`` values.

    Returns:
        Inner join of ``phenotypes`` and ``genotypes`` on normalized ``eid``.

    Raises:
        ValueError: If either frame lacks ``eid``.
        JoinDropRateError: If the unmatched-ID drop rate exceeds ``max_drop_rate``.
    """
    if "eid" not in phenotypes.columns:
        raise ValueError("Phenotype table missing column: eid")
    if "eid" not in genotypes.columns:
        raise ValueError("Genotype matrix missing column: eid")

    pheno_dtype = _eid_dtype_name(phenotypes)
    geno_dtype = _eid_dtype_name(genotypes)
    schema_mismatch = pheno_dtype != geno_dtype
    schema_normalized = schema_mismatch or pheno_dtype != "String" or geno_dtype != "String"

    pheno_norm = _normalize_eid_frame(phenotypes)
    geno_norm = _normalize_eid_frame(genotypes)

    pheno_eids = set(_sorted_unique_eids(pheno_norm))
    geno_eids = set(_sorted_unique_eids(geno_norm))
    matched = pheno_eids & geno_eids
    pheno_only = sorted(pheno_eids - geno_eids)
    geno_only = sorted(geno_eids - pheno_eids)
    union_n = len(pheno_eids | geno_eids)
    drop_rate = (len(pheno_only) + len(geno_only)) / union_n if union_n else 1.0
    halted = drop_rate > max_drop_rate

    if schema_mismatch:
        LOG.warning(
            "eid schema mismatch: phenotype=%s genotype=%s; normalizing both to Utf8",
            pheno_dtype,
            geno_dtype,
        )
    if pheno_only or geno_only:
        LOG.warning(
            "eid join orphans: pheno_only=%d geno_only=%d drop_rate=%.4f",
            len(pheno_only),
            len(geno_only),
            drop_rate,
        )

    if audit_log is not None:
        write_join_audit_log(
            audit_log,
            pheno_n=len(pheno_eids),
            geno_n=len(geno_eids),
            matched_n=len(matched),
            pheno_only=pheno_only,
            geno_only=geno_only,
            pheno_dtype=pheno_dtype,
            geno_dtype=geno_dtype,
            schema_mismatch=schema_mismatch,
            schema_normalized=schema_normalized,
            drop_rate=drop_rate,
            max_drop_rate=max_drop_rate,
            halted=halted,
        )

    if halted:
        raise JoinDropRateError(
            f"eid join drop rate {drop_rate:.4%} exceeds {max_drop_rate:.4%} "
            f"(pheno_only={len(pheno_only)}, geno_only={len(geno_only)}, "
            f"matched={len(matched)}, schema_mismatch={schema_mismatch}; "
            f"phenotype_dtype={pheno_dtype}, genotype_dtype={geno_dtype}). "
            f"Audit log: {audit_log}"
        )

    return pheno_norm.join(geno_norm, on="eid", how="inner")


def ad_diagnosis_from_code(code: str | None) -> int:
    """Binary AD flag: non-empty ICD-style code → 1, else 0 (mock v2 field).

    Args:
        code: Mock AD diagnosis code string, or ``None``.

    Returns:
        ``1`` when ``code`` is non-empty after stripping; otherwise ``0``.
    """
    if code is None:
        return 0
    return 1 if str(code).strip() else 0


def genotype_phenotype_contingency(
    genotype: np.ndarray,
    outcome: np.ndarray,
) -> np.ndarray:
    """Build a 2×3 table: rows = outcome (0, 1), columns = genotype dosage (0, 1, 2).

    Non-finite genotype or outcome values (e.g. missing dosages) are dropped; they
    are never coerced to dosage ``0``.

    Args:
        genotype: Per-sample alt-allele dosages (expected values 0, 1, or 2).
        outcome: Per-sample binary outcomes (expected values 0 or 1).

    Returns:
        ``int64`` array of shape ``(2, 3)`` with counts for valid genotype/outcome pairs.
    """
    g = np.asarray(genotype, dtype=float)
    y = np.asarray(outcome, dtype=float)
    valid = (
        np.isfinite(g)
        & np.isfinite(y)
        & (g >= 0)
        & (g <= 2)
        & (y >= 0)
        & (y <= 1)
    )
    g_i = g[valid].astype(np.int64)
    y_i = y[valid].astype(np.int64)
    if g_i.size == 0:
        return np.zeros((2, 3), dtype=np.int64)
    return np.bincount(y_i * 3 + g_i, minlength=6).reshape(2, 3).astype(np.int64)


def _dominant_2x2_from_contingency(table_2x3: np.ndarray) -> list[list[int]]:
    """Collapse 2×3 genotype table to 2×2 dominant model (carrier = dosage 1 or 2).

    Args:
        table_2x3: Contingency from :func:`genotype_phenotype_contingency`.

    Returns:
        Nested list ``[[control_non_carrier, control_carrier],
        [case_non_carrier, case_carrier]]``.
    """
    control_non_carrier = int(table_2x3[0, 0])
    carrier_control = int(table_2x3[0, 1] + table_2x3[0, 2])
    case_non_carrier = int(table_2x3[1, 0])
    carrier_case = int(table_2x3[1, 1] + table_2x3[1, 2])
    return [
        [control_non_carrier, carrier_control],
        [case_non_carrier, carrier_case],
    ]


def _or_confidence_interval(
    a: float,
    b: float,
    c: float,
    d: float,
    or_point: float,
    *,
    alpha: float = 0.05,
) -> tuple[float, float]:
    """Woolf log-OR interval on a 2×2 table ``[[a,b],[c,d]]``.

    Zero cells receive a Haldane (+0.5) correction for the standard-error
    calculation only. The CI is centered on the caller-supplied ``or_point``
    (typically the same Haldane-adjusted OR returned as the point estimate).

    Args:
        a: Control non-carriers.
        b: Control carriers.
        c: Case non-carriers.
        d: Case carriers.
        or_point: Odds-ratio point estimate used as the CI center.
        alpha: Two-sided type I error rate (default 0.05 → 95% CI).

    Returns:
        ``(ci_low, ci_high)`` on the odds-ratio scale.
    """
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    log_or = math.log(or_point)
    se = math.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)
    z = float(norm.ppf(1.0 - alpha / 2.0))
    return math.exp(log_or - z * se), math.exp(log_or + z * se)


def dominant_odds_ratio(
    genotype: np.ndarray,
    outcome: np.ndarray,
) -> tuple[float, float, float, float, int]:
    """Dominant-model OR (carrier vs non-carrier) with two-sided Fisher exact p and 95% CI.

    Args:
        genotype: Per-sample alt-allele dosages (0 / 1 / 2).
        outcome: Per-sample binary phenotype (0 / 1).

    Returns:
        Tuple ``(OR, CI_low, CI_high, p_value, n)`` where ``n`` is the number of
        samples in the 2×3 contingency table. Degenerate tables may yield
        ``nan`` / ``inf`` odds ratios and undefined CIs.
    """
    table_2x3 = genotype_phenotype_contingency(genotype, outcome)
    n = int(table_2x3.sum())
    table_2x2 = _dominant_2x2_from_contingency(table_2x3)
    a, b = table_2x2[0]
    c, d = table_2x2[1]
    if n == 0:
        return float("nan"), float("nan"), float("nan"), float("nan"), 0

    if a == 0 and b == 0:
        or_val = float("inf") if c > 0 else float("nan")
        p_val = 1.0 if c == 0 else 0.0
        ci_low, ci_high = float("nan"), float("nan")
        return or_val, ci_low, ci_high, p_val, n

    if c == 0 and d == 0:
        or_val = 0.0 if b > 0 else float("nan")
        p_val = 1.0 if b == 0 else 0.0
        ci_low, ci_high = float("nan"), float("nan")
        return or_val, ci_low, ci_high, p_val, n

    _fisher_or, p_val = fisher_exact(table_2x2)
    # When any cell is zero, report the Haldane-adjusted OR so the Woolf CI
    # is centered on the same point estimate that is returned.
    if min(a, b, c, d) == 0:
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        or_float = (aa * dd) / (bb * cc)
    else:
        or_float = float(_fisher_or)
    if not math.isfinite(or_float) or or_float <= 0:
        ci_low, ci_high = float("nan"), float("nan")
    else:
        ci_low, ci_high = _or_confidence_interval(float(a), float(b), float(c), float(d), or_float)
    return or_float, ci_low, ci_high, float(p_val), n


def _snp_columns(frame: pl.DataFrame) -> list[str]:
    return [col for col in frame.columns if col != "eid" and col not in _PHENOTYPE_JOIN_COLS]


_PHENOTYPE_JOIN_COLS: Final[frozenset[str]] = frozenset(
    {
        "eid",
        "age",
        "sex",
        "parental_longevity",
        "ad_diagnosis_code",
        "pd_diagnosis_code",
        "frailty_weight_loss",
        "frailty_exhaustion",
        "frailty_weakness",
        "frailty_slowness",
        "frailty_low_activity",
    }
)


def run_association_scan(
    joined: pl.DataFrame,
    *,
    phenotype_col: PhenotypeColumn,
    outcome: pl.Series | None = None,
) -> pl.DataFrame:
    """Scan each LA-SNP column for dominant-model association with a binary phenotype.

    Args:
        joined: Phenotype–genotype frame with ``eid``, phenotype fields, and SNP columns.
        phenotype_col: Name of the binary phenotype column used for labelling / logging.
            Used as the outcome source when ``outcome`` is ``None``.
        outcome: Optional precomputed binary outcome series (e.g. derived AD flag).
            When set, overrides reading ``phenotype_col`` from ``joined``.

    Returns:
        Association results with columns ``rsID``, ``OR``, ``CI_low``, ``CI_high``,
        ``p_value``, and ``n``.

    Raises:
        ValueError: If ``outcome`` is ``None`` and ``phenotype_col`` is missing, or if
            no SNP columns are present in ``joined``.
    """
    if outcome is None:
        if phenotype_col not in joined.columns:
            raise ValueError(f"Phenotype column missing: {phenotype_col}")
        y = joined[phenotype_col].to_numpy()
    else:
        y = outcome.to_numpy()

    snp_cols = _snp_columns(joined)
    if not snp_cols:
        raise ValueError("No SNP columns found in joined frame")

    rows: list[dict[str, float | str | int]] = []
    for rs_id in sorted(snp_cols):
        g = joined[rs_id].to_numpy()
        contingency = genotype_phenotype_contingency(g, y)
        or_val, ci_low, ci_high, p_val, n = dominant_odds_ratio(g, y)
        LOG.debug(
            "%s × %s 2×3 contingency:\n%s → dominant OR=%.4g p=%.4g n=%d",
            rs_id,
            phenotype_col,
            contingency,
            or_val,
            p_val,
            n,
        )
        rows.append(
            {
                "rsID": rs_id,
                "OR": or_val,
                "CI_low": ci_low,
                "CI_high": ci_high,
                "p_value": p_val,
                "n": n,
            }
        )

    return pl.DataFrame(rows).select(list(LA_SNP_ASSOC_COLUMNS))


def write_association_results(frame: pl.DataFrame, path: Path) -> None:
    """Write association scan CSV with a synthetic-data disclaimer header.

    Args:
        frame: Association results (typically from :func:`run_association_scan`).
        path: Destination CSV path; parent directories are created as needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    body = frame.write_csv()
    path.write_text(f"# {SYNTHETIC_DISCLAIMER}\n{body}", encoding="utf-8")


def run_integration_pipeline(
    pheno_path: Path,
    vcf_path: Path,
    output_dir: Path,
    *,
    audit_log: Path = DEFAULT_AUDIT_LOG,
    max_drop_rate: float = MAX_JOIN_DROP_RATE,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Load mock inputs, join on ``eid``, run LA-SNP scans, write association CSVs under ``output_dir``.

    Args:
        pheno_path: Mock UKB phenotype CSV path.
        vcf_path: Mock LA-SNP VCF path with sample IDs equal to ``eid``.
        output_dir: Directory for ``assoc_la_snp_*.csv`` outputs.
        audit_log: Path for dropped-ID / schema mismatch audit logging.
        max_drop_rate: Maximum allowed unmatched ``eid`` fraction (default 1%).

    Returns:
        Tuple ``(joined, parental_longevity_results, ad_results)``.
    """
    phenotypes = load_phenotype_table(pheno_path)
    genotypes = load_genotype_matrix_from_vcf(vcf_path)
    joined = join_phenotypes_genotypes(
        phenotypes,
        genotypes,
        audit_log=audit_log,
        max_drop_rate=max_drop_rate,
    )

    ad_outcome = pl.Series(
        "ad_diagnosis",
        [ad_diagnosis_from_code(v) for v in joined["ad_diagnosis_code"].to_list()],
    )

    parental_results = run_association_scan(
        joined, phenotype_col="parental_longevity"
    )
    ad_results = run_association_scan(joined, phenotype_col="ad_diagnosis_code", outcome=ad_outcome)

    write_association_results(parental_results, output_dir / PARENTAL_LONGEVITY_OUT)
    write_association_results(ad_results, output_dir / AD_OUT)

    LOG.info(
        "Activity %s: joined %d samples × %d LA-SNPs; wrote %s and %s (%s)",
        ACTIVITY_ID,
        joined.height,
        len(_snp_columns(joined)),
        output_dir / PARENTAL_LONGEVITY_OUT,
        output_dir / AD_OUT,
        SYNTHETIC_DISCLAIMER,
    )
    return joined, parental_results, ad_results
