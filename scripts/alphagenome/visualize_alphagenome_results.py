from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
ALPHAGENOME_DATA_DIR = REPO_ROOT / "analysis" / "alphagenome"
FIGURES_DIR = REPO_ROOT / "figures" / "alphagenome"
IMPACT_CSV = ALPHAGENOME_DATA_DIR / "alphagenome_impact_analysis.csv"
BAR_PLOT = FIGURES_DIR / "alphagenome_impact_bar_plot.png"
SCATTER_PLOT = FIGURES_DIR / "alphagenome_ref_vs_alt_scatter.png"


def create_visualizations():
    try:
        df = pd.read_csv(IMPACT_CSV)
    except FileNotFoundError:
        print(f"Analysis file not found at {IMPACT_CSV}. Run the analysis script first.")
        return

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # Filter for top 20 variants by absolute percentage change
    top_df = df.sort_values(by='abs_perc_change', ascending=False).head(20)

    # 1. Bar plot of percentage changes
    plt.figure(figsize=(12, 8))
    # Create a column for display: "Gene (SNP)"
    top_df['display_name'] = top_df['gene'] + " (" + top_df['snp'] + ")"
    
    # Use a custom color palette based on positive/negative change
    colors = ['#ff7f0e' if x > 0 else '#1f77b4' for x in top_df['perc_change']]
    
    sns.barplot(data=top_df, x='perc_change', y='display_name', palette=colors)
    plt.axvline(x=0, color='black', linestyle='-', linewidth=1)
    plt.title('Top 20 Predicted Regulatory Impacts (RNA-seq % Change)', fontsize=15)
    plt.xlabel('Predicted Expression Change (%)', fontsize=12)
    plt.ylabel('Gene (SNP ID)', fontsize=12)
    plt.tight_layout()
    plt.savefig(BAR_PLOT)
    print(f"Saved {BAR_PLOT}")

    # 2. Scatter plot of Ref vs Alt scores
    plt.figure(figsize=(10, 8))
    plt.scatter(df['ref_score'], df['alt_score'], alpha=0.6, color='purple')
    
    # Add a diagonal line for y=x (no change)
    max_val = max(df['ref_score'].max(), df['alt_score'].max())
    plt.plot([0, max_val], [0, max_val], 'r--', alpha=0.8, label='No Change (y=x)')
    
    plt.title('Predicted RNA-seq Scores: Reference vs Alternate Alleles', fontsize=15)
    plt.xlabel('Reference Allele Score', fontsize=12)
    plt.ylabel('Alternate Allele Score', fontsize=12)
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.7)
    
    # Annotate top variants
    for i, row in top_df.head(5).iterrows():
        plt.annotate(row['gene'], (row['ref_score'], row['alt_score']), 
                     textcoords="offset points", xytext=(0,10), ha='center', fontsize=9)
                     
    plt.tight_layout()
    plt.savefig(SCATTER_PLOT)
    print(f"Saved {SCATTER_PLOT}")

if __name__ == "__main__":
    create_visualizations()
