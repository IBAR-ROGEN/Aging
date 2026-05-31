# ROGEN Aging — Activity map

Maps IBAR-ROGEN activity IDs to code, scripts, and documentation.

| Activity | Description | Package / module | Primary CLI | Documentation |
|----------|-------------|------------------|-------------|---------------|
| **2.1.7.1** | AlphaGenome LA-SNP regulatory comparison | — | `scripts/alphagenome/alphagenome_sequence_comparer.py` | [ALPHAGENOME_ANALYSIS_EXPLANATION.md](ALPHAGENOME_ANALYSIS_EXPLANATION.md) |
| **2.1.7.1** | LA-SNP pathway network figure | — | `scripts/figures/generate_network_fig.py` | [FIGURES.md](FIGURES.md) |
| **2.1.8.1** | Methylation calling pipeline | `rogen_aging.methylation_visualizations` | `pipeline_validation.sh`, `downstream_analysis.R` | [METHYLATION_PIPELINE_README.md](METHYLATION_PIPELINE_README.md) |
| **2.1.8.1** | LA-SNP manifest + public AF validation | `rogen_aging.ukb` | `rogen-ukb-manifest`, `rogen-compare-af-gnomad` | [LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md) |
| **2.1.8.1** | Synthetic UKB-RAP mock folder | `rogen_aging.vcf`, `scripts/ukb/mock_rap_folder.py` | `scripts/ukb/mock_rap_folder.py` | [SYNTHETIC_UKB_RAP_GENERATOR.md](SYNTHETIC_UKB_RAP_GENERATOR.md) |
| **2.1.11.1** | Synthetic UKB integrative validation | `rogen_aging.integration` | `rogen-ukb-integrate` | [UKB_INTEGRATION_PIPELINE.md](UKB_INTEGRATION_PIPELINE.md) |
| — | Epigenetic clock library | `rogen_aging.clock` | `rogen-clock` | [CLOCK_LIBRARY.md](CLOCK_LIBRARY.md) |
| — | Multi-omics EDA dashboard | `rogen_aging.eda_dashboard` | Streamlit app | [EDA_DASHBOARD.md](EDA_DASHBOARD.md) |

## Script layout (by workflow)

```
scripts/
├── clock/           run_clock.py, validate_clock.py (deprecated wrappers), Romanian demo
├── ukb/             manifest, gnomad compare, mock generators, integration
├── vcf/             synthetic Romanian cohort VCF
├── figures/         matplotlib / networkx manuscript renders
├── alphagenome/     sequence comparer + analysis
├── eda/             mock epigenetic EDA
└── dev/             security hook, CI audit, R bootstrap, utilities
```

Deprecated one-line shims at the old flat `scripts/*.py` paths forward to these folders.
