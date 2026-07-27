# Nomenclature reconciliation & manuscript figures

**Script:** [`reconcile_and_generate_figures.py`](../reconcile_and_generate_figures.py)  
**Related:** [FIGURES.md](FIGURES.md) · [AF_COMPARISON_FIGURES.md](AF_COMPARISON_FIGURES.md) · [GENOMICS_ANALYSIS.md](GENOMICS_ANALYSIS.md)

## Purpose

End-to-end manuscript finalization helper that:

1. **Audits gene nomenclature** — maps historical CETP / HLA / HSP70 locus labels in manuscript text and tables to HGNC symbols, then cross-checks the 41-gene candidate list against Supplementary Table 3 (symbols, variant counts, functional clusters).
2. **Plots allele frequencies** — log-scale scatter of `ROGEN_AF` vs gnomAD v4 NFE (`gnomAD_AF`), highlighting loci with `|ΔAF| > 0.05`.
3. **Draws the 41-gene network** — STRING edges laid out by functional pathway; nodes coloured by Pro-/Anti-/Context-Dependent longevity class and sized by degree centrality.

## Required inputs

| Path | Role |
|------|------|
| `manuscript/tables/Supplementary_Table_3.xlsx` | 41-gene summary (`Gene_Symbol`, `Variant_Count`, `Functional_Cluster`, optional `Longevity_Class`) |
| `manuscript/tables/41_gene_candidate_list.csv` | Candidate list used for concordance |
| `manuscript/text/*.md` | Draft text scanned for legacy locus names |
| `analysis/rogen_vs_gnomad_af.csv` | Columns `rsID`, `ROGEN_AF`, `gnomAD_AF` (aliases `AF_1kg` / `AF_gnomad_nfe` accepted) |
| `data/network/41_gene_interactions.csv` | STRING edges (`gene_a`, `gene_b`, `string_score`) |
| `data/network/41_gene_nodes.csv` | Optional node metadata (falls back to Supplementary Table 3) |

## Outputs

| Path | Description |
|------|-------------|
| `outputs/nomenclature_audit.log` | Legacy → HGNC reconciliations and concordance findings |
| `outputs/figures/Figure_AF_Scatter.pdf` / `.png` | Log-scale AF scatter (300 dpi PNG + vector PDF) |
| `outputs/figures/Figure_41_Gene_Network.pdf` / `.png` | Pathway-grouped STRING network |

Exit code `1` if any concordance **DISCREPANCY** remains after the audit.

## CLI usage

```bash
uv sync

uv run python reconcile_and_generate_figures.py
```

Override paths when needed:

```bash
uv run python reconcile_and_generate_figures.py \
  --supp-table manuscript/tables/Supplementary_Table_3.xlsx \
  --candidate-list manuscript/tables/41_gene_candidate_list.csv \
  --text-dir manuscript/text \
  --af-csv analysis/rogen_vs_gnomad_af.csv \
  --network-csv data/network/41_gene_interactions.csv \
  --fig-dir outputs/figures \
  --audit-log outputs/nomenclature_audit.log
```

```bash
uv run python reconcile_and_generate_figures.py --help
```
