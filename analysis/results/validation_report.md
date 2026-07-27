# Genomics validation report (GRCh38/hg38)

Independent recomputation from `overlapping_genes_with_snps.xlsx`. Manuscript-stated values are not trusted; both stated and computed values are recorded.

- API cache hits: 66
- API cache misses: 0

## Stated vs computed

| Metric | Stated | Computed | Verdict | Notes |
|--------|--------|----------|---------|-------|
| Unique gene symbols (SNP Association == significant) [GRCh38/hg38] | 41 | 41 | MATCH |  |
| Unique gene symbols (SNP Association == non-significant) [GRCh38/hg38] | NA | 53 | NEEDS REVIEW | Genes in both groups (11): ASIC2, CETP, GPX1, HLA-DQB1, HSPA1A, HSPA1B, HSPA1L, LRP1B, RAD23B, SDC4, TLR4 |
| LA-SNP count (a) significant raw rows [GRCh38/hg38] | 70 | 79 | MISMATCH |  |
| LA-SNP count (b) unique (gene, SNP Identifier) pairs [GRCh38/hg38] | 70 | 64 | MISMATCH |  |
| LA-SNP count (c) unique (gene, canonical rsID) excluding HLA/repeat [GRCh38/hg38] | 70 | 48 | MISMATCH | Unique canonical rsIDs across genes: 48 |
| Unique canonical rsIDs (significant, non-HLA/repeat) [GRCh38/hg38] | 58 | 48 | MISMATCH | Manuscript/AlphaGenome docs cite 58 unique rsIDs for the LA-SNP set |
| AlphaGenome LA-SNP rows (external curated table) [GRCh38/hg38] | 70 | 70 | MATCH | Source: analysis/alphagenome/alphagenome_impact_analysis.csv |
| AlphaGenome unique rsIDs [GRCh38/hg38] | 58 | 58 | MATCH |  |
| CETP unique SNP Identifier count (significant, raw names) [GRCh38/hg38] | 10 | 10 | MATCH | identifiers: CETP, I405V, rs1800774, rs1800777, rs289714, rs4784744, rs5882, rs5883, rs9923854, rs9930761 |
| CETP unique canonical rsID count (significant) [GRCh38/hg38] | 10 | 9 | MISMATCH | canonical rsIDs: rs1273184461, rs1800774, rs1800777, rs289714, rs4784744, rs5882, rs5883, rs9923854, rs9930761 |
| HSPA1A unique SNP Identifier count (significant, raw names) [GRCh38/hg38] | 3 | 3 | MATCH | identifiers: -110A>C, G190C, rs1008438 |
| HSPA1A unique canonical rsID count (significant) [GRCh38/hg38] | 3 | 1 | MISMATCH | canonical rsIDs: rs1008438 |
| HSPA1B unique SNP Identifier count (significant, raw names) [GRCh38/hg38] | 3 | 3 | MATCH | identifiers: 1267A/G, 1267A>G, A1267G |
| HSPA1B unique canonical rsID count (significant) [GRCh38/hg38] | 3 | 0 | MISMATCH | canonical rsIDs: NA |
| HSPA1L unique SNP Identifier count (significant, raw names) [GRCh38/hg38] | 2 | 2 | MATCH | identifiers: 2437T>C, T2437C |
| HSPA1L unique canonical rsID count (significant) [GRCh38/hg38] | 2 | 0 | MISMATCH | canonical rsIDs: NA |
| NDUFS1 unique SNP Identifier count (significant, raw names) [GRCh38/hg38] | 5 | 5 | MATCH | identifiers: rs11548670, rs11695633, rs1801318, rs6435324, rs6435326 |
| NDUFS1 unique canonical rsID count (significant) [GRCh38/hg38] | 5 | 5 | MATCH | canonical rsIDs: rs11548670, rs11695633, rs1801318, rs6435324, rs6435326 |
| PCSK1 unique SNP Identifier count (significant, raw names) [GRCh38/hg38] | 3 | 3 | MATCH | identifiers: rs155979, rs3762985, rs3811952 |
| PCSK1 unique canonical rsID count (significant) [GRCh38/hg38] | 3 | 3 | MATCH | canonical rsIDs: rs155979, rs3762985, rs3811952 |
| APOC1 unique SNP Identifier count (significant, raw names) [GRCh38/hg38] | 2 | 2 | MATCH | identifiers: APOC1, rs4420638 |
| APOC1 unique canonical rsID count (significant) [GRCh38/hg38] | 2 | 1 | MISMATCH | canonical rsIDs: rs4420638 |
| HLA-DQB1 unique SNP Identifier count (significant, raw names) [GRCh38/hg38] | 2 | 2 | MATCH | identifiers: DQB103, DQB105 |
| HLA-DQB1 unique canonical rsID count (significant) [GRCh38/hg38] | 2 | 0 | MISMATCH | canonical rsIDs: NA |
| SDC4 unique SNP Identifier count (significant, raw names) [GRCh38/hg38] | 2 | 2 | MATCH | identifiers: rs1981429, rs2251252 |
| SDC4 unique canonical rsID count (significant) [GRCh38/hg38] | 2 | 2 | MATCH | canonical rsIDs: rs1981429, rs2251252 |
| Duplicate alias check CETP I405V vs rs5882 [GRCh38/hg38] | same variant (manuscript convention) | different rsIDs | MISMATCH | CETP: I405V -> ['rs1273184461']; rs5882 -> ['rs5882']; same_variant=False on GRCh38/hg38 (Ensembl-resolved rsIDs differ) |
| Duplicate alias check HSPA1A -110A>C vs rs1008438 [GRCh38/hg38] | same variant (manuscript convention) | different rsIDs | NEEDS REVIEW | HSPA1A: -110A>C -> ['NA']; rs1008438 -> ['rs1008438']; unresolved |
| APOC1 row labeled 'APOC1' [GRCh38/hg38] | real SNP | gene-name placeholder (HpaI RFLP study row) | NEEDS REVIEW | Row uses gene symbol as SNP Identifier for an RFLP study; canonical rsID=NA unless resolved from literature/dbSNP |
| Cluster list size AD-up [GRCh38/hg38] | 410 | NA | NEEDS REVIEW | Supplementary Table 3.xlsx not in repository |
| Cluster list size AD-down [GRCh38/hg38] | 833 | NA | NEEDS REVIEW | Supplementary Table 3.xlsx not in repository |
| Cluster list size PD-up [GRCh38/hg38] | 318 | NA | NEEDS REVIEW | Supplementary Table 3.xlsx not in repository |
| Cluster list size PD-down [GRCh38/hg38] | 229 | NA | NEEDS REVIEW | Supplementary Table 3.xlsx not in repository |
| All 41 significant genes present in >=1 cluster list [GRCh38/hg38] | yes | NA | NEEDS REVIEW | Missing file: /Users/mityatoren/code/rogen_aging/data/Supplementary Table 3.xlsx |

## Proposed corrections

### LA-SNP count 70

Manuscript cites 70 pairs / 58 unique rsIDs (AlphaGenome table matches). The overlap xlsx has 79 significant rows and 64 unique (gene, identifier) pairs on GRCh38/hg38 — inflated by duplicate rsIDs (e.g. rs4420638×3), legacy aliases, and gene-name placeholders.

### CETP I405V vs rs5882

On GRCh38/hg38, Ensembl maps I405V to rs1273184461 (chr16:56983397) and rs5882 to chr16:56982180 — different loci. Manuscript likely treats them as synonymous; table double-counts if both are kept.

### HSPA1A -110A>C

Ensembl HGVS recoder cannot resolve c.-110A>C to an rsID; rs1008438 is listed separately at chr6:31815431. Proposed merge pending dbSNP synonym confirmation (currently NA for -110A>C).

### Gene placeholders

Rows using gene symbols as SNP Identifier (APOC1, CETP, PRR5L, SGK1, VEGFA, YWHAG) are not dbSNP rsIDs — replace with specific variants or exclude from SNP counts.

### HLA / HMOX1

DQB103, DQB105, and (GT)n repeat are not SNPs — exclude from LA-SNP totals (already excluded in definition c).

### Supplementary Table 3

Add data/Supplementary Table 3.xlsx to verify all 41 significant genes appear in AD/PD cluster lists.


## Flagged issues

- **gene_coordinates**: 124 rows have Start > End on GRCh38/hg38; using min/max span for checks
- **non_snp_exclusion**: Excluded from LA-SNP SNP counts: (GT)n repeat
- **non_snp_exclusion**: Excluded from LA-SNP SNP counts: DQB103
- **non_snp_exclusion**: Excluded from LA-SNP SNP counts: DQB105
- **duplicate_alias**: CETP: I405V -> ['rs1273184461']; rs5882 -> ['rs5882']; same_variant=False on GRCh38/hg38 (Ensembl-resolved rsIDs differ)
- **duplicate_alias**: HSPA1A: -110A>C -> ['NA']; rs1008438 -> ['rs1008438']; unresolved
- **gene_placeholder**: APOC1 significant row uses 'APOC1' as SNP Identifier — likely not a dbSNP rsID
- **coordinate**: APOC1/rs4420638: GRCh38 pos 44919689 outside span 44914324-44919345 (GRCh38/hg38)
- **coordinate**: APOC1/rs4420638: GRCh38 pos 44919689 outside span 44914324-44919345 (GRCh38/hg38)
- **coordinate**: APOC1/rs4420638: GRCh38 pos 44919689 outside span 44914324-44919345 (GRCh38/hg38)
- **coordinate**: HSPA1A/rs1008438: GRCh38 pos 31815431 outside span 31815542-31817941 (GRCh38/hg38)
- **coordinate**: HSPA1A/rs1008438: GRCh38 pos 31815431 outside span 31815542-31817941 (GRCh38/hg38)
- **coordinate**: HSPA1A/rs1008438: GRCh38 pos 31815431 outside span 31815542-31817941 (GRCh38/hg38)
- **coordinate**: MT2A/rs28366003: GRCh38 pos 56608579 outside span 56608583-56609496 (GRCh38/hg38)
- **coordinate**: PCSK1/rs155979: GRCh38 pos 96434194 outside span 96390332-96433247 (GRCh38/hg38)
- **coordinate**: PCSK1/rs3762985: GRCh38 pos 96435023 outside span 96390332-96433247 (GRCh38/hg38)
- **coordinate**: PCSK1/rs155979: GRCh38 pos 96434194 outside span 96390332-96433247 (GRCh38/hg38)
- **coordinate**: PCSK1/rs3762985: GRCh38 pos 96435023 outside span 96390332-96433247 (GRCh38/hg38)
- **coordinate**: SEMA4A/rs1468772: GRCh38 pos 156146697 outside span 156147372-156177743 (GRCh38/hg38)
- **cluster_lists**: Supplementary Table 3 not found at /Users/mityatoren/code/rogen_aging/data/Supplementary Table 3.xlsx; cannot verify 41 significant genes against AD/PD cluster lists (GRCh38/hg38)
