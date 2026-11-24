# Local eQTL Catalogue tables — provenance

Verification only. These tables were built outside this repo
(see [INTEGRATIVE_PIPELINE.md](INTEGRATIVE_PIPELINE.md): there is no
in-repo writer for `results/eqtl_associations_full.csv` /
`results/annotation_layer.csv`). This document records what is known
about the local copies and the result of running
`scripts/eqtl/verify_local_eqtl_tables.py`.

## What is known

- Source: tables were built outside this repo (see
  `docs/INTEGRATIVE_PIPELINE.md`).
- File mtime on `results/eqtl_associations_full.csv` and
  `results/annotation_layer.csv` is **2026-08-30** and is **not
  authoritative** (filesystem timestamp only).
- Locked dataset list: **14** `dataset_id` values in
  `analysis/eqtl/eqtl_datasets_used.tsv`.
- Catalogue release version: **unknown, not recoverable from the repo**.

## Dataset IDs

- `QTD000146`
- `QTD000151`
- `QTD000156`
- `QTD000161`
- `QTD000166`
- `QTD000171`
- `QTD000176`
- `QTD000181`
- `QTD000186`
- `QTD000191`
- `QTD000196`
- `QTD000201`
- `QTD000206`
- `QTD000356`

## Inventory (from verification)

| Path | sha256 | n_rows | mtime (UTC) | rows on 47 panel rsIDs | panel rsIDs present |
|------|--------|--------|-------------|------------------------|---------------------|
| `data/processed/variants_47_input_v2.csv` | `c3d1393a18943d2f99227f90186ed88a1ef896a5df451b77cc0d6432ec7be4e5` | 47 | 2026-10-29T10:50:25Z | 47 | 47 |
| `analysis/eqtl/eqtl_datasets_used.tsv` | `cbdd03ea30c80b80006a199564107e4e845c4fd87c4b32a6ad391bd687296cb7` | 14 | 2026-10-28T19:16:53Z | n/a | n/a |
| `results/eqtl_associations_full.csv` | `acdb93d21949e146db4c39113a10ae0ab951ceeea8f0f03e40bb27291ed66b90` | 1650 | 2026-08-30T17:34:00Z | 1650 | 43 |
| `results/annotation_layer.csv` | `d2aad734cf1d4005bd1db0d545b6a4dfb2f7b3607c0242bf6fd2ca22ad138afd` | 73 | 2026-08-30T17:34:00Z | 73 | 47 |

### Columns

- `data/processed/variants_47_input_v2.csv`: `chrom`, `pos`, `ref`, `alt`, `rsid`, `gene_symbol`, `alt_source`, `alt_evidence`, `alt_status`
- `analysis/eqtl/eqtl_datasets_used.tsv`: `dataset_id`
- `results/eqtl_associations_full.csv`: `canonical_rsid`, `eqtl_alt_allele`, `tissue`, `dataset_id`, `n_samples`, `gene_id`, `egene_symbol`, `is_cis_to_annotated_gene`, `beta`, `se`, `pvalue`, `fdr_bh`, `maf`, `median_tpm`
- `results/annotation_layer.csv`: `gene`, `canonical_rsid`, `chrom`, `pos_grch38`, `ref_allele`, `alt_allele`, `n_alleles`, `is_multiallelic`, `category`, `most_severe_consequence`, `AF_1000G_EUR`, `AF_1000G_EUR_status`, `AF_gnomADe_NFE`, `AF_gnomADe_NFE_status`, `AF_gnomADg_NFE`, `AF_gnomADg_NFE_status`, `AF_gnomAD_v4_NFE`, `AF_gnomAD_v4_NFE_status`, `delta_AF`, `vep_consequence_allele`, `vep_impact`, `n_transcripts`, `canonical_transcript`, `canonical_consequence`, `vep_gene_symbols`, `n_regulatory_features`, `regulatory_biotypes`, `n_eqtl_tests`, `n_tissues_tested`, `n_eqtl_fdr05`, `n_tissues_fdr05`, `tissues_fdr05`, `top_egene`, `top_egene_id`, `top_tissue`, `top_beta`, `top_se`, `top_pvalue`, `top_fdr_bh`, `top_is_cis`, `any_cis_fdr05`, `cis_egenes_fdr05`, `is_ensembl_minor`, `ensembl_global_maf`, `alphamissense_eligible`, `risk_allele`, `risk_allele_source`, `eqtl_query_status`

## Integrated evidence check (temporary output)

`scripts/longevity/integrated_evidence.py` was invoked with the local
associations, annotation, and dataset list, writing under a temporary
directory (committed `analysis/longevity/` outputs were not written).

- Variants in temp table: **47**
- Labels: GWAS **1** / eQTL **16** / population **5**
- Ground truth compared: GWAS **1** / eQTL **16** / population **5**

## Verification provenance

- Git commit: `430f0558a0ad899c0f7dab0ff63fc7128ee0bf71`
- Timestamp (UTC): `2026-11-24T09:22:39Z`
- Python: `3.12.8`
- Polars: `1.35.1`
- Typer: `0.24.1`
- Machine-readable sidecar: `analysis/eqtl/verify_local_eqtl_tables_provenance.json`
