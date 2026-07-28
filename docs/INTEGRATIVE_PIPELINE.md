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
| Annotated variants | `chrom`, `pos`, `ref`, `alt` | Production default: July `Combined_Master` (GTEx summary cols stripped, then rebuilt) |
| Long eQTL table | `rsid`, `tissue`, `nes`, `p_value` | Production default: July `GTEx_eQTL_Summary` (`slope`→`nes`); aliases `rsID` accepted |
| AlphaGenome scores (optional) | `rsid` or coordinates | Usually already present on Combined_Master |
| Probe annotation (optional) | `IlmnID`, `UCSC_RefGene_Name` | Semicolon-separated gene lists exploded on join |
| Sample genotypes (optional) | `sample_id`, `rsid`, `alt_dosage` | Extra phenotype columns are carried through aggregation |

Typical upstream sources: July annotation workbook
(`outputs/Supplementary_Table_1_Annotated_Variants.xlsx`) or parquet siblings
`Supplementary_Table_1_Combined_Master.parquet` /
`Supplementary_Table_1_GTEx_eQTL_Summary.parquet`.

## Outputs

Written under `analysis/integrative/results/` by the production CLIs (Parquet):

| Artefact | Contents |
|----------|----------|
| `annotated_variants.parquet` | Variants + GTEx summary (+ optional AlphaGenome) |
| `eqtl_summary.parquet` | One row per rsID: best tissue, n hits, tissue list |
| `variant_risks.parquet` | Channel scores + `composite_risk` ∈ [0, 1] |
| `methylation_links.parquet` | Variant–probe links via gene symbol (optional) |
| `sample_profiles.parquet` | Dosage-weighted `sample_risk` per sample (optional) |

Demo / fixture smoke tests write under `analysis/integrative/demo/` when `--demo` is passed.

## CLI usage

```bash
uv sync --extra dev

# Production (Activity A.2.1.11.1): July Combined_Master × GTEx → results/
# Requires outputs/Supplementary_Table_1_Annotated_Variants.xlsx (or parquet siblings).
uv run python scripts/integrative/run_pipeline.py
# → analysis/integrative/results/{annotated_variants,eqtl_summary,variant_risks}.parquet

# Tissue map only (same production defaults)
uv run python scripts/integrative/map_variant_tissues.py

# Risk scores from the production tissue map
uv run python scripts/integrative/integrate_phenotypes.py

# Offline fixture smoke test
uv run python scripts/integrative/run_pipeline.py --demo
```

## Library usage

```python
import polars as pl
from rogen_aging.integrative import run_integrative_pipeline
from rogen_aging.integrative.io import load_production_eqtls, load_production_variants

result = run_integrative_pipeline(
    load_production_variants(),
    load_production_eqtls(),
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
