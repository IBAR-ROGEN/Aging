# AlphaGenome Analysis Explanation

This document provides a detailed explanation of the methodology and logic behind the "AlphaGenome Sequence Comparer" script and its regulatory impact predictions.

## 1. Overview of the Prediction
The analysis identifies how a single genetic change (a SNP) is predicted to "re-program" the expression of a specific gene. The primary output is a **predicted regulatory score**, which represents the **abundance of RNA transcripts** (gene expression) in a specific tissue.

For this analysis, we targeted the **Brain (Caudate nucleus, UBERON:0001157)** tissue.

## 2. Methodology: How it Works
AlphaGenome is a deep learning model trained on the human genome and various functional genomics datasets (like CAGE, DNase-seq, and RNA-seq). It has "learned" the language of DNA and can predict how a given sequence will be transcribed into RNA by analyzing its regulatory landscape (enhancers, promoters, etc.).

### Step-by-Step Process:
1.  **Coordinate Mapping**: We use genomic coordinates (GRCh38) for specific genes of interest.
2.  **Sequence Fetching**: The script fetches a fixed **128kb (131,072 bp)** window of DNA sequence from the Ensembl REST API, centered on the target SNP position.
3.  **Variant Retrieval**: The script fetches the Reference and Alternate (Longevity-Associated) alleles for the target SNP from the Ensembl Variation API.

## 3. The Comparison Logic (Reference vs. Alternate)

### The "Control" (Reference)
*   **Sequence**: The standard human reference genome sequence (GRCh38).
*   **Action**: The model predicts the expression level based on this "natural" or wild-type sequence context.

### The "Comparison" (Alternate)
*   **Sequence**: A "mutated" version of the reference sequence.
*   **Action**: The script swaps the single nucleotide at the SNP position from the Reference allele to the **Longevity-Associated Alternate allele**.
*   **Prediction**: The model re-predicts the expression level for this virtually mutated sequence.

## 4. Interpreting the Results
The final result is expressed as a **Percentage Change (%)** between the Reference and Alternate predictions. This represents the model's estimate of the **functional impact** of the longevity variant.

*   **Positive % Change**: The alternate allele is predicted to **increase** (upregulate) the expression of that gene in the target tissue.
*   **Negative % Change**: The alternate allele is predicted to **decrease** (downregulate/knock down) the expression of that gene.

### Biological Example: NDUFS1
In the mitochondrial gene **NDUFS1**, the model predicted that the longevity-associated SNP **rs6435324** reduces the gene's expression in the brain by approximately **4.78%**. This suggests that a slight reduction in this mitochondrial component might be a mechanism linked to increased longevity.

## 5. Technical Details
*   **Script**: `scripts/alphagenome_sequence_comparer.py`
*   **Analysis**: `scripts/analyze_alphagenome_results.py`
*   **Input Data**: `overlapping_genes_with_snps.xlsx`
*   **Results File**: `alphagenome_comparison_results.csv`
*   **Processed Impact**: `alphagenome_impact_analysis.csv`
*   **Target Length**: 131,072 bp (128kb)
*   **Ontology Term**: UBERON:0001157 (Caudate nucleus)
