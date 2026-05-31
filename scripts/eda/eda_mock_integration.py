"""
Exploratory Data Analysis (EDA) for mock clinical epigenetic aging data.

This script loads a hypothetical CSV with clinical and epigenetic measurements,
computes Epigenetic Age Acceleration (EAA) residuals, prints descriptive stats,
and generates a scatter plot of chronological vs epigenetic age.

Expected CSV columns: Sample_ID, Chronological_Age, Epigenetic_Age, Phenotype_Score
"""

from pathlib import Path

import pandas as pd
import seaborn as sns

# ------------------------------------------------------------------------------
# Configuration: paths relative to project root
# ------------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
INPUT_CSV = PROJECT_ROOT / "test_data" / "mock_epigenetic_clinical.csv"
OUTPUT_DIR = PROJECT_ROOT / "results"
OUTPUT_PLOT = OUTPUT_DIR / "mock_eaa_plot.png"

# ------------------------------------------------------------------------------
# 1. Load data and compute EAA residuals
# ------------------------------------------------------------------------------
# Epigenetic Age Acceleration (EAA) = Epigenetic_Age - Chronological_Age
# Positive values indicate "accelerated aging"; negative values indicate "slower" aging.
df = pd.read_csv(INPUT_CSV)
df["EAA_Residuals"] = df["Epigenetic_Age"] - df["Chronological_Age"]

# ------------------------------------------------------------------------------
# 2. Print basic descriptive statistics
# ------------------------------------------------------------------------------
print("=" * 60)
print("Descriptive Statistics (mock clinical data)")
print("=" * 60)
print(df.describe())
print()
print("Row count:", len(df))
print("Missing values per column:")
print(df.isnull().sum())
print("=" * 60)

# ------------------------------------------------------------------------------
# 3. Scatter plot: Chronological Age vs Epigenetic Age
# ------------------------------------------------------------------------------
# Ensure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Create scatter plot with a diagonal reference line (y = x)
# Points above the line = accelerated aging; below = decelerated.
sns.set_style("whitegrid")
ax = sns.scatterplot(
    data=df,
    x="Chronological_Age",
    y="Epigenetic_Age",
    alpha=0.7,
    edgecolor="none",
)
# Add diagonal reference line (perfect agreement between clocks)
ax.axline((0, 0), slope=1, color="gray", linestyle="--", alpha=0.7, label="y = x")
ax.set_title("Chronological Age vs Epigenetic Age")
ax.set_xlabel("Chronological Age")
ax.set_ylabel("Epigenetic Age")
ax.legend()

ax.figure.tight_layout()
ax.figure.savefig(OUTPUT_PLOT, dpi=150, bbox_inches="tight")
print(f"Plot saved to: {OUTPUT_PLOT}")
