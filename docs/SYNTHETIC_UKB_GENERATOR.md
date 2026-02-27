# Synthetic UK Biobank Data Generator

**Project:** IBAR-ROGEN Aging  
**Script:** `scripts/mock_ukb_generator.py`  
**Output:** `test_data/mock_clinical_data.csv` (default)

## Overview

The synthetic UK Biobank data generator produces fake tabular data that mimics UK Biobank format. It is intended for pipeline development and testing **before real data arrives**, and is safe to commit to GitHub — no real participant identifiers or sensitive data are used.

## Columns

| Column | Description | Values |
|--------|-------------|--------|
| `Sample_ID` | Synthetic sample identifier | `MOCK_00000001`, `MOCK_00000002`, … |
| `Age` | Age in years | 40–80 (uniform) |
| `EAA` | Epigenetic Age Acceleration | Normal(0, 5) |
| `rs_mock_001` … `rs_mock_005` | Dummy SNP genotypes | 0 (ref/ref), 1 (ref/alt), 2 (alt/alt) |

SNP genotypes follow a Hardy–Weinberg distribution with minor allele frequency (MAF) 0.2 by default.

## Usage

### Basic (1000 samples, default output)

```bash
uv run scripts/mock_ukb_generator.py
```

Writes to `test_data/mock_clinical_data.csv`.

### Custom sample count and output path

```bash
uv run scripts/mock_ukb_generator.py --n-samples 500 --output data/synthetic_cohort.csv
```

### Reproducibility

Use `--seed` for reproducible random data:

```bash
uv run scripts/mock_ukb_generator.py --seed 42
```

Use `--seed 0` to disable the seed (non-deterministic).

### All options

```bash
uv run scripts/mock_ukb_generator.py --help
```

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--n-samples` | `-n` | 1000 | Number of synthetic samples |
| `--output` | `-o` | `test_data/mock_clinical_data.csv` | Output CSV path |
| `--seed` | `-s` | 42 | Random seed (0 = no seed) |

## Example Output

```csv
Sample_ID,Age,EAA,rs_mock_001,rs_mock_002,rs_mock_003,rs_mock_004,rs_mock_005
MOCK_00000001,43,7.17,0,0,1,0,0
MOCK_00000002,71,0.46,1,0,0,1,0
MOCK_00000003,66,2.90,1,0,0,1,1
...
```

## Security and Compliance

- **Sample_ID prefix:** `MOCK_` (not `UKB_`) so the pre-commit hook does not block commits.
- **Whitelisted:** The script and `test_data/mock_clinical_data.csv` are excluded from the UK Biobank security content scan.
- **No real data:** All values are synthetic. Safe for public repositories.

## Use Cases

1. **Pipeline development** — Test analysis code before Romanian or UK Biobank data is available.
2. **CI/CD** — Run integration tests on synthetic data without real participant data.
3. **Documentation** — Demonstrate expected input formats for downstream tools.
4. **Tutorials** — Provide reproducible example data for training.

## Related Documentation

- [UK Biobank Pre-Commit Hook](UKB_PRE_COMMIT_HOOK.md) — Security checks and whitelisting
- [Project Structure](PROJECT_STRUCTURE.md) — Where `test_data/` and `data/` fit in the layout

---

**Last Updated:** February 27, 2026
