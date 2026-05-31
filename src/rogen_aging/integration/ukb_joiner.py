"""Synthetic UKB phenotype–genotype join and LA-SNP association scan.

Activity **2.1.11.1** — integrative architecture validation on **synthetic** mock RAP
output only (from ``scripts/ukb_mock_gen.py``). No real UK Biobank data; outputs are
for pipeline QA, **not** biological conclusions.
"""

from __future__ import annotations

import logging
import math
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
    """Load mock UKB phenotype CSV (``eid`` + v2 fields) with safety comment lines skipped."""
    if not path.is_file():
        raise FileNotFoundError(f"Phenotype table not found: {path.resolve()}")
    return pl.read_csv(path, comment_prefix="#")


def _alt_dosage_from_gt_type(gt_type: int) -> int | None:
    """Map cyvcf2 ``gt_type`` to alt-allele dosage 0/1/2; unknown genotypes → ``None``."""
    if gt_type == 0:
        return 0
    if gt_type == 1:
        return 1
    if gt_type == 2:
        return 2
    return None


def load_genotype_matrix_from_vcf(path: Path) -> pl.DataFrame:
    """Load LA-SNP VCF into a wide matrix: one row per ``eid``, one column per rsID (0/1/2)."""
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
        dosages: list[int | None] = []
        for gt_type in variant.gt_types:
            dosages.append(_alt_dosage_from_gt_type(int(gt_type)))
        dosage_by_snp[str(rs_id)] = dosages

    if not dosage_by_snp:
        raise ValueError(f"VCF contains no variant records: {path}")

    frame = pl.DataFrame({"eid": sample_ids})
    for rs_id, values in dosage_by_snp.items():
        frame = frame.with_columns(pl.Series(rs_id, values))
    return frame


def join_phenotypes_genotypes(
    phenotypes: pl.DataFrame,
    genotypes: pl.DataFrame,
) -> pl.DataFrame:
    """Inner-join phenotype table and genotype matrix on ``eid`` (one row per participant)."""
    if "eid" not in phenotypes.columns:
        raise ValueError("Phenotype table missing column: eid")
    if "eid" not in genotypes.columns:
        raise ValueError("Genotype matrix missing column: eid")

    joined = phenotypes.join(genotypes, on="eid", how="inner")
    if joined.height != phenotypes.height:
        raise ValueError(
            f"Join row count {joined.height} != phenotype rows {phenotypes.height}; "
            "eid sets may not match."
        )
    return joined


def ad_diagnosis_from_code(code: str | None) -> int:
    """Binary AD flag: non-empty ICD-style code → 1, else 0 (mock v2 field)."""
    if code is None:
        return 0
    return 1 if str(code).strip() else 0


def genotype_phenotype_contingency(
    genotype: np.ndarray,
    outcome: np.ndarray,
) -> np.ndarray:
    """Build a 2×3 table: rows = outcome (0, 1), columns = genotype dosage (0, 1, 2)."""
    table = np.zeros((2, 3), dtype=np.int64)
    g = np.asarray(genotype, dtype=np.int64)
    y = np.asarray(outcome, dtype=np.int64)
    valid = (g >= 0) & (g <= 2) & (y >= 0) & (y <= 1)
    g = g[valid]
    y = y[valid]
    for outcome_val in (0, 1):
        for dosage in (0, 1, 2):
            table[outcome_val, dosage] = int(np.sum((y == outcome_val) & (g == dosage)))
    return table


def _dominant_2x2_from_contingency(table_2x3: np.ndarray) -> list[list[int]]:
    """Collapse 2×3 genotype table to 2×2 dominant model (carrier = dosage 1 or 2)."""
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
    """Woolf log-OR interval on a 2×2 table [[a,b],[c,d]] with Haldane correction if needed."""
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        or_point = (a * d) / (b * c)
    log_or = math.log(or_point)
    se = math.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)
    z = float(norm.ppf(1.0 - alpha / 2.0))
    return math.exp(log_or - z * se), math.exp(log_or + z * se)


def dominant_odds_ratio(
    genotype: np.ndarray,
    outcome: np.ndarray,
) -> tuple[float, float, float, float, int]:
    """Dominant-model OR (carrier vs non-carrier) with two-sided Fisher exact p and 95% CI."""
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

    or_val, p_val = fisher_exact(table_2x2)
    or_float = float(or_val)
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
    """Scan each LA-SNP column for dominant-model association with a binary phenotype."""
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
    """Write association scan CSV with a synthetic-data disclaimer header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    body = frame.write_csv()
    path.write_text(f"# {SYNTHETIC_DISCLAIMER}\n{body}", encoding="utf-8")


def run_integration_pipeline(
    pheno_path: Path,
    vcf_path: Path,
    output_dir: Path,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Load mock inputs, join on ``eid``, run LA-SNP scans, write association CSVs under ``output_dir``."""
    phenotypes = load_phenotype_table(pheno_path)
    genotypes = load_genotype_matrix_from_vcf(vcf_path)
    joined = join_phenotypes_genotypes(phenotypes, genotypes)

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
