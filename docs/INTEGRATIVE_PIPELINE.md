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

STALE: the production July workbook and parquet siblings still use alleles from `data/processed/variants_47_input.csv`. They were not rebuilt after `variants_47_input_v2.csv` because the July stage joins AlphaGenome. Affected rsids: rs28366003, rs1800774, rs41383, rs9916344, rs7207422, rs1801318, rs1981429, rs139170, rs155979, rs524533, and rs882696.

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

## Integrated public evidence (per-variant labels)

After the longevity GWAS lookup (Prompt 4) has written
`analysis/longevity/longevity_gwas_lookup.csv`, join that table with the
September eQTL Catalogue layer, the panel frequency/TSI–NFE test, and the
illustrative power grid. Join keys are **rsid + panel alt**.

```bash
uv run python scripts/longevity/integrated_evidence.py
```

| Input | Role |
|-------|------|
| `analysis/longevity/longevity_gwas_lookup.csv` | GWAS direction and BH q per trait on the panel alt |
| `results/eqtl_associations_full.csv` | September eQTL Catalogue associations (alt allele, `dataset_id`, `beta`, `fdr_bh`) |
| `analysis/eqtl/eqtl_datasets_used.tsv` | Locked list of 14 `dataset_id` values |
| `analysis/panel/ro_longevity_panel.csv` | Gene, frequencies, TSI vs NFE q / `tsi_differs_fdr05` |
| `analysis/panel/ro_longevity_power_grid.csv` | Illustrative power grid for the power tier |

| Output | Contents |
|--------|----------|
| `analysis/longevity/integrated_evidence.csv` | One row per panel variant |
| `analysis/longevity/integrated_evidence_dot_matrix.png` | Dot matrix: variants × evidence labels, sorted by number of labels |
| `analysis/longevity/integrated_evidence_provenance.json` | Input checksums, label definitions, release-check record |
| `analysis/longevity/release_checks/integrated_evidence.json` | `rogen-release-check` summary report |

### Evidence labels (definitions)

A variant can carry several labels. **No label implies causality**, Romanian
enrichment, or a ROGEN association.

| Label | Definition |
|-------|------------|
| `GWAS support` | BH q < 0.05 in at least one longevity/lifespan trait in the GWAS lookup |
| `eQTL support` | FDR < 0.05 in at least one eQTL Catalogue dataset for the panel alt |
| `population difference` | TSI vs NFE BH q < 0.05 on the panel alt (public reference panels) |

`eqtl_alt_effect_direction` is `consistent` when significant datasets agree on
the sign of the mean alt-allele beta, `mixed` when signs disagree, and `none`
when there are no significant datasets.

`power_tier` uses the illustrative rank cell (1000 cases per arm, allelic OR
1.5, Bonferroni \(0.05/47\)): **high** ≥ 0.8, **moderate** ≥ 0.5, **low** <
0.5, **unavailable** when power is missing.

## Tests

```bash
uv run pytest tests/test_integrative.py tests/test_integrated_evidence.py -q
```

Covers chromosome/rsID normalisation, multi-allelic `variant_key`, eQTL
filtering/summaries, AlphaGenome joins, methylation links, channel scores,
composite-risk bounds, sample aggregation, CLI-facing pipeline wiring, and the
integrated evidence label rules.
