# Integrative multi-omics pipeline (variants × tissues × phenotype)

**Package:** `rogen_aging.integrative`  
**Scripts:** `scripts/integrative/`  
**Related:** [JULY_ANNOTATION_PIPELINE.md](JULY_ANNOTATION_PIPELINE.md) · [LA_SNP_GTEX_ANNOTATION.md](LA_SNP_GTEX_ANNOTATION.md) · [LA_SNP_VEP_ANNOTATION.md](LA_SNP_VEP_ANNOTATION.md) · [UKB_INTEGRATION_PIPELINE.md](UKB_INTEGRATION_PIPELINE.md)

## Overview

Offline joining and scoring layer that cross-references molecularly annotated
sequencing variants with tissue-specific GTEx eQTL profiles, optional
methylation probe links, and composite phenotypic risk. Live GTEx / VEP API
fetching stays in the annotation scripts; this module only consumes tables.

Core API:

| Symbol | Role |
|--------|------|
| `VariantTissueMapper` | Normalise loci, summarise eQTLs, join AlphaGenome, link probes |
| `PhenotypeIntegrator` | Channel scores → `composite_risk` → optional sample profiles |
| `run_integrative_pipeline` | End-to-end mapper + integrator convenience wrapper |

Default GTEx tissues are brain regions plus whole blood
(`DEFAULT_TARGET_TISSUES`). Risk channel weights default to VEP 0.25,
AlphaGenome 0.25, AlphaMissense 0.25, GTEx 0.15, epigenetic 0.10
(`DEFAULT_WEIGHTS`).

## Inputs

| Table | Required columns | Notes |
|-------|------------------|-------|
| Annotated variants | `chrom`, `pos`, `ref`, `alt` | Optional: `rsid`, `gene_symbol`, `vep_impact`, Alpha*/GTEx columns, `age_acceleration` |
| Long eQTL table | `rsid`, `tissue`, `nes`, `p_value` | Optional: `gene_symbol` / `eqtl_gene_symbol`, `gtex_variant_id` |
| AlphaGenome scores (optional) | `rsid` or coordinates | Aliases `diff`, `perc_change`, `snp`, … normalised to `alphagenome_*` |
| Probe annotation (optional) | `IlmnID`, `UCSC_RefGene_Name` | Semicolon-separated gene lists exploded on join |
| Sample genotypes (optional) | `sample_id`, `rsid`, `alt_dosage` | Extra phenotype columns are carried through aggregation |

Typical upstream sources: July annotation workbook / GTEx CSV under
`analysis/gtex_annotation/`, AlphaGenome matrices under `data/scores/` or
`analysis/alphagenome/`.

## Outputs

Written under `analysis/integrative/` by the CLIs (Parquet):

| Artefact | Contents |
|----------|----------|
| `annotated` / `annotated_variants.parquet` | Variants + GTEx summary (+ optional AlphaGenome) |
| `eqtl_summary.parquet` | One row per rsID: best tissue, n hits, tissue list |
| `variant_risks.parquet` | Channel scores + `composite_risk` ∈ [0, 1] |
| `methylation_links.parquet` | Variant–probe links via gene symbol (optional) |
| `sample_profiles.parquet` | Dosage-weighted `sample_risk` per sample (optional) |

## CLI usage

```bash
uv sync --extra dev

# End-to-end: tissue map + composite risk (+ optional samples / probes)
uv run python scripts/integrative/run_pipeline.py \
  --variants data/processed/prioritized_variants.csv \
  --eqtls analysis/gtex_annotation/la_snp_gtex_eqtls.csv \
  --alphagenome data/scores/alphagenome_raw.parquet \
  --probes data/annotation/hm450_probe_genes.csv \
  --samples data/processed/sample_genotypes.csv \
  --output-dir analysis/integrative/

# Tissue map only
uv run python scripts/integrative/map_variant_tissues.py \
  --variants data/processed/prioritized_variants.csv \
  --eqtls analysis/gtex_annotation/la_snp_gtex_eqtls.csv \
  -o analysis/integrative/

# Risk scores from an already-mapped table
uv run python scripts/integrative/integrate_phenotypes.py \
  --annotated analysis/integrative/annotated_variants.parquet \
  --samples data/processed/sample_genotypes.csv \
  -o analysis/integrative/
```

## Library usage

```python
import polars as pl
from rogen_aging.integrative import run_integrative_pipeline

result = run_integrative_pipeline(
    pl.read_csv("data/processed/prioritized_variants.csv"),
    pl.read_csv("analysis/gtex_annotation/la_snp_gtex_eqtls.csv"),
)
risks = result["variant_risks"].sort("composite_risk", descending=True)
```

## Tests

```bash
uv run pytest tests/test_integrative.py -q
```

Covers chromosome/rsID normalisation, multi-allelic `variant_key`, eQTL
filtering/summaries, AlphaGenome joins, methylation links, channel scores,
composite-risk bounds, sample aggregation, and CLI-facing pipeline wiring.
