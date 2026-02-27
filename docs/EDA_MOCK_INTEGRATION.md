# EDA Mock Integration Script

**Project:** IBAR-ROGEN Aging  
**Script:** `scripts/eda_mock_integration.py`  
**Input:** `test_data/mock_epigenetic_clinical.csv`  
**Output:** `results/mock_eaa_plot.png`

## Overview

The EDA mock integration script performs exploratory data analysis on hypothetical clinical epigenetic aging data. It loads a CSV with chronological age, epigenetic age, and phenotype scores; computes Epigenetic Age Acceleration (EAA) residuals; prints descriptive statistics; and generates a scatter plot.

This script is intended for:
- Testing analysis pipelines before real data is available
- Demonstrating expected input/output formats for downstream tools
- Integration testing of the epigenetic aging workflow

## Expected CSV Columns

| Column              | Description                          |
|---------------------|--------------------------------------|
| `Sample_ID`         | Sample or subject identifier         |
| `Chronological_Age` | Age in years (calendar age)          |
| `Epigenetic_Age`    | Age predicted by methylation clock   |
| `Phenotype_Score`   | Optional phenotype/outcome measure   |

## Outputs

1. **EAA_Residuals** — Computed as `Epigenetic_Age - Chronological_Age`:
   - Positive = accelerated epigenetic aging
   - Negative = decelerated epigenetic aging

2. **Descriptive statistics** — Printed to stdout (count, mean, std, min, quartiles, max, missing values)

3. **Scatter plot** — Saved to `results/mock_eaa_plot.png`:
   - X-axis: Chronological Age
   - Y-axis: Epigenetic Age
   - Diagonal reference line (y = x) for perfect agreement

## Usage

```bash
uv run python scripts/eda_mock_integration.py
```

Ensure `test_data/mock_epigenetic_clinical.csv` exists. A sample file is versioned in the repository.

## Example Input

```csv
Sample_ID,Chronological_Age,Epigenetic_Age,Phenotype_Score
S001,45,48.2,0.72
S002,52,50.1,0.55
S003,38,42.3,0.81
...
```

## Dependencies

- pandas
- seaborn
- matplotlib (via seaborn)

All are listed in `pyproject.toml`.
