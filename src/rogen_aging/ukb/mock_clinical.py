"""Synthetic UK Biobank-style tabular clinical data generator."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import typer

DUMMY_SNP_IDS: tuple[str, ...] = (
    "rs_mock_001",
    "rs_mock_002",
    "rs_mock_003",
    "rs_mock_004",
    "rs_mock_005",
)


def generate_synthetic_ukb_data(
    n_samples: int = 1000,
    min_age: int = 40,
    max_age: int = 80,
    eaa_mean: float = 0.0,
    eaa_std: float = 5.0,
    snp_maf: float = 0.2,
    seed: int | None = 42,
) -> pd.DataFrame:
    """Generate synthetic UK Biobank-style tabular data."""
    rng = np.random.default_rng(seed)
    sample_ids = [f"MOCK_{i:08d}" for i in range(1, n_samples + 1)]
    age = rng.integers(min_age, max_age + 1, size=n_samples)
    sex = rng.integers(0, 2, size=n_samples)
    bmi = rng.uniform(15.0, 50.0, size=n_samples)
    ad_diagnosis = rng.binomial(1, 0.02, size=n_samples)
    eaa = rng.normal(loc=eaa_mean, scale=eaa_std, size=n_samples)

    p_0 = (1 - snp_maf) ** 2
    p_1 = 2 * snp_maf * (1 - snp_maf)
    p_2 = snp_maf**2
    probs = [p_0, p_1, p_2]

    snp_cols = {
        snp_id: rng.choice([0, 1, 2], size=n_samples, p=probs) for snp_id in DUMMY_SNP_IDS
    }

    return pd.DataFrame(
        {
            "Sample_ID": sample_ids,
            "Age": age,
            "Sex": sex,
            "BMI": bmi,
            "AD_diagnosis": ad_diagnosis,
            "EAA": eaa,
            **snp_cols,
        }
    )


app = typer.Typer(
    help="Generate synthetic UK Biobank-style mock clinical data for pipeline testing."
)


@app.command()
def main(
    n_samples: int = typer.Option(
        1000,
        "--n-samples",
        "-n",
        help="Number of synthetic samples to generate.",
    ),
    output: Path = typer.Option(
        Path("test_data/mock_clinical_data.csv"),
        "--output",
        "-o",
        path_type=Path,
        help="Output CSV path.",
    ),
    seed: int | None = typer.Option(
        42,
        "--seed",
        "-s",
        help="Random seed for reproducibility. Use 0 for no seed.",
    ),
) -> None:
    """Generate synthetic UK Biobank-style mock clinical data."""
    resolved_seed = None if seed == 0 else seed
    df = generate_synthetic_ukb_data(n_samples=n_samples, seed=resolved_seed)
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    typer.echo(f"Wrote {len(df)} synthetic samples to {output}")
