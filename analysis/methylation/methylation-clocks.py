"""
Horvath 2013 epigenetic clock implementation.

Predicts DNA methylation age (DNAm age) from Illumina array beta values
(450K, 27K, or EPIC) using the 353-CpG multi-tissue clock.

==== !!!! IMPORTANT !!!! ==== we need to adapt this to the ROGEN dataset, which deals with Nanopore data, not Illumina arrays

Reference:
  Horvath S. (2013). DNA methylation age of human tissues and cell types.
  Genome Biology, 14(10), R115. https://doi.org/10.1186/gb-2013-14-10-r115

Coefficients: Use the official file from
  https://horvath.genetics.ucla.edu/html/dnamage/
  or Table S2 from the paper supplementary material.
"""

from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd


# Default path for coefficient file (same directory as this module)
# horvath_353_coefficients.csv for an alternate title
# currently, gb-2013-14-10-r115-S3.csv is the file that contains the coefficients for the Horvath 2013 epigenetic clock (Additional file 3 from the paper
# Supplementary Table S3)
_DEFAULT_COEF_PATH = Path(__file__).resolve().parent / "test_data" / "gb-2013-14-10-r115-S3.csv"


def _load_horvath_coefficients(path: Union[str, Path, None] = None) -> tuple[float, dict[str, float]]:
    """
    Load intercept and per-probe coefficients from a CSV.

    CSV format: one column 'Probe' (or 'probe_id') and one 'Coefficient' (or 'coef').
    First row can be header. Optionally a row with Probe='Intercept' for intercept;
    otherwise intercept is assumed 0.

    Returns:
        (intercept, dict of probe_id -> coefficient)
    """
    path = path or _DEFAULT_COEF_PATH
    path = Path(path)
    intercept = 0.0
    coefs = {}

    if not path.exists():
        raise FileNotFoundError(
            f"Horvath coefficient file not found: {path}\n"
            "Download the 353 CpG coefficients from https://horvath.genetics.ucla.edu/html/dnamage/"
            " or from the paper supplementary (Table S2), and save as CSV with columns 'Probe' and 'Coefficient'."
        )

    # gb-2013-14-10-r115-S3.csv has 2 leading comment/empty lines, then CpGmarker,CoefficientTraining,...
    with open(path, encoding="utf-8") as f:
        first_line = f.readline()
    skip = 2 if "CpGmarker" not in first_line and "probe" not in first_line.lower() else 0
    df = pd.read_csv(path, skiprows=range(skip) if skip else None)
    probe_col = next(
        (c for c in df.columns if "probe" in c.lower() or "cg" in c.lower() or c == "CpGmarker"),
        df.columns[0],
    )
    coef_col = next(
        (c for c in df.columns if "coefficient" in c.lower() or "coef" in c.lower() or "weight" in c.lower()),
        df.columns[1],
    )

    for _, row in df.iterrows():
        probe = str(row[probe_col]).strip()
        if not probe or probe.startswith("#"):
            continue
        try:
            val = float(row[coef_col])
        except (TypeError, ValueError):
            continue
        if probe.lower().strip("()") == "intercept":
            intercept = val
        else:
            coefs[probe] = val

    return intercept, coefs


class HorvathClock:
    """
    Horvath 2013 epigenetic clock: 353 CpG sites, multi-tissue.

    Input: beta values (0–1) keyed by Illumina probe ID (e.g. cg16867657).
    Missing CpGs are imputed with 0.5 (standard choice).
    """

    def __init__(self, coefficient_path: Union[str, Path, None] = None):
        self.intercept, self.coefficients = _load_horvath_coefficients(coefficient_path)
        self.probe_ids = list(self.coefficients.keys())

    def _single_sample_age(self, beta: dict[str, float], missing_impute: float = 0.5) -> float:
        """Predict DNAm age for one sample given probe_id -> beta."""
        raw = self.intercept
        for probe, coef in self.coefficients.items():
            b = beta.get(probe)
            if b is None or (isinstance(b, float) and np.isnan(b)):
                b = missing_impute
            raw += coef * float(b)
        # Clip to plausible human age range (Horvath output can exceed calendar age)
        return float(np.clip(raw, 0.0, 120.0))

    def predict(
        self,
        beta: Union[dict[str, float], pd.DataFrame],
        missing_impute: float = 0.5,
    ) -> Union[float, np.ndarray]:
        """
        Predict DNAm age (Horvath clock).

        Parameters
        ----------
        beta : dict or DataFrame
            - dict: single sample, probe_id -> beta value (0–1).
            - DataFrame: rows = probes (index or column named 'Probe'/'probe_id'),
              columns = samples; or columns = probes, rows = samples (one row).
        missing_impute : float
            Value to use for missing CpGs (default 0.5).

        Returns
        -------
        float or ndarray
            Predicted DNAm age in years (one value or one per sample).
        """
        if isinstance(beta, dict):
            return self._single_sample_age(beta, missing_impute)

        df = pd.DataFrame(beta)
        # DataFrame: assume either (probes x samples) or (samples x probes)
        if df.index.astype(str).str.match(r"cg\d+", na=False).any():
            # Rows = probes
            probe_axis, sample_axis = 0, 1
        elif df.columns.astype(str).str.match(r"cg\d+", na=False).any():
            # Columns = probes
            probe_axis, sample_axis = 1, 0
        else:
            raise ValueError(
                "DataFrame must have probe IDs (e.g. cg16867657) as index or column names."
            )

        if probe_axis == 0:
            # rows = probes
            ages = []
            for col in df.columns:
                ages.append(self._single_sample_age(df[col].to_dict(), missing_impute))
            return np.array(ages)
        else:
            # columns = probes
            ages = []
            for _, row in df.iterrows():
                ages.append(self._single_sample_age(row.to_dict(), missing_impute))
            return np.array(ages)

    def __len__(self) -> int:
        return len(self.coefficients)


def horvath_clock(
    beta: Union[dict[str, float], pd.DataFrame],
    coefficient_path: Union[str, Path, None] = None,
    missing_impute: float = 0.5,
) -> Union[float, np.ndarray]:
    """
    Convenience function: predict DNAm age using the Horvath 2013 clock.

    Parameters
    ----------
    beta : dict or DataFrame
        Methylation beta values keyed by probe ID (see HorvathClock.predict).
    coefficient_path : path or None
        Path to CSV with columns Probe, Coefficient. If None, uses
        horvath_353_coefficients.csv next to this file.
    missing_impute : float
        Value for missing CpGs (default 0.5).

    Returns
    -------
    float or ndarray
        Predicted DNAm age in years.
    """
    clock = HorvathClock(coefficient_path=coefficient_path)
    return clock.predict(beta, missing_impute=missing_impute)


# Path to Aging/test_data (repo test_data folder)
_TEST_DATA_DIR = Path(__file__).resolve().parent.parent.parent / "test_data"
_COEF_FILE = "gb-2013-14-10-r115-S3.csv"


def _load_test_betas_from_csv(path: Path, beta_col: str = "medianByCpG") -> dict[str, float]:
    """Load methylation beta values from the coefficient CSV (e.g. medianByCpG, medianByCpGYoung, medianByCpGOld)."""
    df = pd.read_csv(path, skiprows=range(2))
    probe_col = "CpGmarker"
    if beta_col not in df.columns:
        raise ValueError(f"Column {beta_col!r} not in CSV. Available: {list(df.columns)}")
    beta = {}
    for _, row in df.iterrows():
        probe = str(row[probe_col]).strip()
        if not probe or probe.lower().strip("()") == "intercept":
            continue
        val = row[beta_col]
        if pd.isna(val):
            continue
        try:
            beta[probe] = float(val)
        except (TypeError, ValueError):
            continue
    return beta


if __name__ == "__main__":
    coef_path = _TEST_DATA_DIR / _COEF_FILE
    if not coef_path.exists():
        print(f"ERROR: Coefficient file not found: {coef_path}")
        print("Expected: Aging/test_data/gb-2013-14-10-r115-S3.csv")
        raise SystemExit(1)

    print("Horvath 2013 clock — test using data from test_data/gb-2013-14-10-r115-S3.csv")
    print("-" * 50)
    print(f"Coefficient file: {coef_path}")

    clock = HorvathClock(coefficient_path=coef_path)
    print(f"Loaded {len(clock)} CpG coefficients, intercept = {clock.intercept:.4f}")

    # Test: beta values from CSV (medianByCpG = median methylation in training set)
    beta_median = _load_test_betas_from_csv(coef_path, "medianByCpG")
    print(f"Loaded {len(beta_median)} beta values from CSV (medianByCpG)")
    age_median = clock.predict(beta_median, missing_impute=0.5)
    print(f"DNAm age (medianByCpG): {age_median:.2f} years")

    # Test: medianByCpGYoung (median methylation in young subset)
    beta_young = _load_test_betas_from_csv(coef_path, "medianByCpGYoung")
    print(f"Loaded {len(beta_young)} beta values from CSV (medianByCpGYoung)")
    age_young = clock.predict(beta_young, missing_impute=0.5)
    print(f"DNAm age (medianByCpGYoung): {age_young:.2f} years")

    # Test: medianByCpGOld (median methylation in old subset)
    beta_old = _load_test_betas_from_csv(coef_path, "medianByCpGOld")
    print(f"Loaded {len(beta_old)} beta values from CSV (medianByCpGOld)")
    age_old = clock.predict(beta_old, missing_impute=0.5)
    print(f"DNAm age (medianByCpGOld): {age_old:.2f} years")

    print("-" * 50)
    print("Test completed successfully.")
