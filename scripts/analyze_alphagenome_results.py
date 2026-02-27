import pandas as pd
import numpy as np
import ast
import re

def parse_track_data(data_str):
    """
    Parses the string representation of AlphaGenome TrackData to extract mean values.
    """
    try:
        parts = data_str.split('alternate=')
        if len(parts) < 2:
            return None, None
            
        ref_part = parts[0]
        alt_part = parts[1]
        
        # Regex to find nonzero_mean in the table footer
        ref_match = re.search(r'nonzero_mean.*?0\s+([\d\.e\-]+)', ref_part, re.DOTALL)
        alt_match = re.search(r'nonzero_mean.*?0\s+([\d\.e\-]+)', alt_part, re.DOTALL)
        
        ref_score = float(ref_match.group(1)) if ref_match else 0.0
        alt_score = float(alt_match.group(1)) if alt_match else 0.0
        
        return ref_score, alt_score
    except Exception as e:
        return 0.0, 0.0

def analyze_results():
    df = pd.read_csv('alphagenome_comparison_results.csv')
    
    # Filter for successful runs
    success_df = df[df['status'] == 'success'].copy()
    
    if success_df.empty:
        print("No successful AlphaGenome results found in the CSV.")
        return

    # DEBUG: Print the first successful row's output to see the exact format
    # print("DEBUG: First row output snippet:")
    # print(str(success_df.iloc[0]['outputs'])[:1000])

    analysis = []
    for idx, row in success_df.iterrows():
        out_str = str(row['outputs'])
        ref_score, alt_score = parse_track_data(out_str)
        
        # If the scores are 0, let's try to see if they are in scientific notation elsewhere
        if ref_score == 0 and alt_score == 0:
            # Maybe it's not nonzero_mean, let's look for anything that looks like a score
            # The TrackData values array starts with some numbers
            val_match = re.search(r'values=array\(\[\[([\d\.e\-]+)', out_str)
            if val_match:
                # Use the first value as a proxy if we can't find the mean
                ref_score = float(val_match.group(1))
                # For alternate, it's after 'alternate='
                alt_val_match = re.search(r'alternate=.*?values=array\(\[\[([\d\.e\-]+)', out_str)
                alt_score = float(alt_val_match.group(1)) if alt_val_match else ref_score

        if ref_score is not None:
            diff = alt_score - ref_score
            perc_change = (diff / ref_score * 100) if ref_score != 0 else 0
            analysis.append({
                'gene': row['gene'],
                'snp': row['snp'],
                'ref_score': ref_score,
                'alt_score': alt_score,
                'diff': diff,
                'perc_change': perc_change
            })
    
    analysis_df = pd.DataFrame(analysis)
    
    if analysis_df.empty:
        print("Could not extract scores from any successful results. The CSV string representation might be missing metadata.")
        return
    
    # Sort by absolute percentage change to find "meaningful" results
    analysis_df['abs_perc_change'] = analysis_df['perc_change'].abs()
    top_results = analysis_df.sort_values(by='abs_perc_change', ascending=False)
    
    print("Top 10 variants by predicted regulatory impact (RNA-seq change):")
    print(top_results[['gene', 'snp', 'ref_score', 'alt_score', 'perc_change']].head(10).to_string(index=False))
    
    # Save the processed analysis
    top_results.to_csv('alphagenome_impact_analysis.csv', index=False)
    print("\nFull impact analysis saved to alphagenome_impact_analysis.csv")

if __name__ == "__main__":
    analyze_results()
