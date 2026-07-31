"""Synthetic UK Biobank-style tabular clinical data generator."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import typer

from rogen_aging.config import cfg_path, get_config
from rogen_aging.config.cli import config_option, load_cli_config

DUMMY_SNP_IDS: tuple[str, ...] = (
    "rs_mock_001",
    "rs_mock_002",
    "rs_mock_003",
    "rs_mock_004",
    "rs_mock_005",
)


def generate_synthetic_ukb_data(
    n_samples: int | None = None,
    min_age: int | None = None,
    max_age: int | None = None,
    eaa_mean: float = 0.0,
    eaa_std: float | None = None,
    snp_maf: float | None = None,
    seed: int | None = 42,
) -> pd.DataFrame:
    """Generate a synthetic UK Biobank-style clinical table with mock SNP genotypes.

    Args:
        n_samples: Number of synthetic participants. Defaults to ``ukb.mock_n_samples``.
        min_age: Inclusive lower bound for simulated age (years).
        max_age: Inclusive upper bound for simulated age (years).
        eaa_mean: Mean of the simulated epigenetic age acceleration (EAA).
        eaa_std: Standard deviation of simulated EAA.
        snp_maf: Minor allele frequency used for Hardy–Weinberg genotype draws.
        seed: Random seed; ``None`` leaves the generator unseeded.

    Returns:
        DataFrame with ``Sample_ID``, demographics, ``AD_diagnosis``, ``EAA``,
        and columns for each ID in ``DUMMY_SNP_IDS`` (dosages 0/1/2).
    """
    cfg = get_config().ukb
    resolved_n = int(cfg.mock_n_samples if n_samples is None else n_samples)
    resolved_min_age = int(cfg.mock_age_min if min_age is None else min_age)
    resolved_max_age = int(cfg.mock_age_max if max_age is None else max_age)
    resolved_eaa_std = float(cfg.mock_eaa_std if eaa_std is None else eaa_std)
    resolved_maf = float(cfg.mock_maf if snp_maf is None else snp_maf)

    rng = np.random.default_rng(seed)
    sample_ids = [f"MOCK_{i:08d}" for i in range(1, resolved_n + 1)]
    age = rng.integers(resolved_min_age, resolved_max_age + 1, size=resolved_n)
    sex = rng.integers(0, 2, size=resolved_n)
    bmi = rng.uniform(15.0, 50.0, size=resolved_n)
    ad_diagnosis = rng.binomial(1, 0.02, size=resolved_n)
    eaa = rng.normal(loc=eaa_mean, scale=resolved_eaa_std, size=resolved_n)

    p_0 = (1 - resolved_maf) ** 2
    p_1 = 2 * resolved_maf * (1 - resolved_maf)
    p_2 = resolved_maf**2
    probs = [p_0, p_1, p_2]

    snp_cols = {snp_id: rng.choice([0, 1, 2], size=resolved_n, p=probs) for snp_id in DUMMY_SNP_IDS}

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
    config: Path | None = config_option(),
    n_samples: int | None = typer.Option(
        None,
        "--n-samples",
        "-n",
        help="Number of synthetic samples to generate. Default: from config.",
    ),
    output: Path | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Output CSV path. Default: from config.",
    ),
    seed: int | None = typer.Option(
        42,
        "--seed",
        "-s",
        help="Random seed for reproducibility. Use 0 for no seed.",
    ),
) -> None:
    """Generate synthetic UK Biobank-style mock clinical data."""
    load_cli_config(config)
    cfg = get_config()
    resolved_n = int(cfg.ukb.mock_n_samples if n_samples is None else n_samples)
    resolved_output = output or cfg_path(cfg, "paths", "ukb", "mock_clinical")
    resolved_seed = None if seed == 0 else seed
    df = generate_synthetic_ukb_data(n_samples=resolved_n, seed=resolved_seed)
    resolved_output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(resolved_output, index=False)
    typer.echo(f"Wrote {len(df)} synthetic samples to {resolved_output}")
