import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def create_visualizations():
    try:
        df = pd.read_csv('alphagenome_impact_analysis.csv')
    except FileNotFoundError:
        print("Analysis file not found. Run the analysis script first.")
        return

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
    plt.savefig('alphagenome_impact_bar_plot.png')
    print("Saved alphagenome_impact_bar_plot.png")

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
    plt.savefig('alphagenome_ref_vs_alt_scatter.png')
    print("Saved alphagenome_ref_vs_alt_scatter.png")

if __name__ == "__main__":
    create_visualizations()
