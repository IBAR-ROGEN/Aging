# Toolchain inventory — September 2026

Checked again on 2026-09-27 on this machine: macOS 27.0, Apple Silicon (`arm64`).

The CLI tools below are not on the default `PATH`. They are installed in a micromamba environment named `rogen-tools` (micromamba 2.9.0). Invoke them as:

```bash
micromamba run -n rogen-tools samtools --version
```

Do not hardcode the environment's install prefix. Homebrew is still absent.

Python packages were checked with `importlib.metadata` in the project environment (`uv run python`, `.venv`, Python 3.12.8) and in system `python3`. Both environments lack `pysam` and `pod5`.

| Tool | Found | Version | Install command (this OS) |
|------|-------|---------|---------------------------|
| samtools | yes | 1.24 | `micromamba create -y -n rogen-tools -c conda-forge -c bioconda samtools=1.24` |
| bcftools | yes | 1.24 | `micromamba create -y -n rogen-tools -c conda-forge -c bioconda bcftools=1.24` |
| minimap2 | yes | 2.31-r1302 | `micromamba create -y -n rogen-tools -c conda-forge -c bioconda minimap2=2.31` |
| modkit | yes | 0.6.4 | `micromamba create -y -n rogen-tools -c conda-forge -c bioconda ont-modkit=0.6.4` |
| tabix | yes | htslib 1.24 | `micromamba create -y -n rogen-tools -c conda-forge -c bioconda htslib=1.24` |
| bgzip | yes | htslib 1.24 | `micromamba create -y -n rogen-tools -c conda-forge -c bioconda htslib=1.24` |
| pysam | no | — | `uv add pysam` |
| pod5 | no | — | `uv add pod5` |

Versions were read with `micromamba run -n rogen-tools <tool> --version`. `tabix` and `bgzip` come from the `htslib=1.24` package. `modkit` is the Bioconda package `ont-modkit=0.6.4`. One create command installed the set together:

```bash
micromamba create -y -n rogen-tools -c conda-forge -c bioconda \
  samtools=1.24 bcftools=1.24 htslib=1.24 minimap2=2.31 ont-modkit=0.6.4
```
