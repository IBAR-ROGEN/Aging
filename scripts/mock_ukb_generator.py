"""Generate synthetic UK Biobank-style tabular data for pipeline development and testing.

Creates fake data mimicking UK Biobank format: Sample_ID, Age, Epigenetic Age Acceleration
(EAA) scores, and dummy SNP genotype columns (0, 1, 2). Safe for use on GitHub without
real participant data.

Usage:
    uv run scripts/mock_ukb_generator.py
    uv run scripts/mock_ukb_generator.py --n-samples 500 --output test_data/mock_clinical_data.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd
import typer

app = typer.Typer(
    help="Generate synthetic UK Biobank-style mock clinical data for pipeline testing."
)

# Dummy SNP IDs (generic rs-style placeholders for longevity/aging-related genes)
DUMMY_SNP_IDS = ["rs_mock_001", "rs_mock_002", "rs_mock_003", "rs_mock_004", "rs_mock_005"]


def generate_synthetic_ukb_data(
    n_samples: int = 1000,
    min_age: int = 40,
    max_age: int = 80,
    eaa_mean: float = 0.0,
    eaa_std: float = 5.0,
    snp_maf: float = 0.2,
    seed: int | None = 42,
) -> pd.DataFrame:
    """Generate synthetic UK Biobank-style tabular data.

    Args:
        n_samples: Number of synthetic samples (rows).
        min_age: Minimum age in years.
        max_age: Maximum age in years.
        eaa_mean: Mean of Epigenetic Age Acceleration (EAA) distribution.
        eaa_std: Standard deviation of EAA.
        snp_maf: Minor allele frequency for SNP genotypes (0–0.5).
        seed: Random seed for reproducibility.

    Returns:
        DataFrame with Sample_ID, Age, EAA, and 5 dummy SNP columns.
    """
    rng = np.random.default_rng(seed)

    # Sample_ID: synthetic IDs (MOCK_ prefix to satisfy UKB security hooks; mimics UKB format)
    sample_ids = [f"MOCK_{i:08d}" for i in range(1, n_samples + 1)]

    # Age: uniform over [min_age, max_age]
    age = rng.integers(min_age, max_age + 1, size=n_samples)

    # Epigenetic Age Acceleration: normal distribution
    eaa = rng.normal(loc=eaa_mean, scale=eaa_std, size=n_samples)

    # SNP genotypes: 0 (ref/ref), 1 (ref/alt), 2 (alt/alt)
    # Using Hardy–Weinberg approximation: P(0)=(1-maf)^2, P(1)=2*maf*(1-maf), P(2)=maf^2
    p_0 = (1 - snp_maf) ** 2
    p_1 = 2 * snp_maf * (1 - snp_maf)
    p_2 = snp_maf**2
    probs = [p_0, p_1, p_2]

    snp_cols = {}
    for snp_id in DUMMY_SNP_IDS:
        snp_cols[snp_id] = rng.choice([0, 1, 2], size=n_samples, p=probs)

    df = pd.DataFrame(
        {
            "Sample_ID": sample_ids,
            "Age": age,
            "EAA": eaa,
            **snp_cols,
        }
    )

    return df


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
    if seed == 0:
        seed = None

    df = generate_synthetic_ukb_data(n_samples=n_samples, seed=seed)

    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    typer.echo(f"Wrote {len(df)} synthetic samples to {output}")


if __name__ == "__main__":
    app()
