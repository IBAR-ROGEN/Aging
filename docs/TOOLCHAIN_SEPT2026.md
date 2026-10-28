# Toolchain inventory — September 2026

Checked again on 2026-09-27 on this machine: macOS 27.0, Apple Silicon (`arm64`).

The CLI tools below are not on the default `PATH`. They are installed in a micromamba environment named `rogen-tools` (micromamba 2.9.0). Invoke them as:

```bash
micromamba run -n rogen-tools samtools --version
```

Do not hardcode the environment's install prefix. Homebrew is still absent.

CLI tools were checked on 2026-09-27. `pysam` and `pod5` were absent then; they are now an optional project extra `seq` in `pyproject.toml` (`uv sync --extra seq`). Intake still reports `pod5 metadata not scanned` when the `pod5` package is not installed.

| Tool | Found | Version | Install command (this OS) |
|------|-------|---------|---------------------------|
| samtools | yes | 1.24 | `micromamba create -y -n rogen-tools -c conda-forge -c bioconda samtools=1.24` |
| bcftools | yes | 1.24 | `micromamba create -y -n rogen-tools -c conda-forge -c bioconda bcftools=1.24` |
| minimap2 | yes | 2.31-r1302 | `micromamba create -y -n rogen-tools -c conda-forge -c bioconda minimap2=2.31` |
| modkit | yes | 0.6.4 | `micromamba create -y -n rogen-tools -c conda-forge -c bioconda ont-modkit=0.6.4` |
| tabix | yes | htslib 1.24 | `micromamba create -y -n rogen-tools -c conda-forge -c bioconda htslib=1.24` |
| bgzip | yes | htslib 1.24 | `micromamba create -y -n rogen-tools -c conda-forge -c bioconda htslib=1.24` |
| pysam | optional (`seq`) | `>=0.22` | `uv sync --extra seq` |
| pod5 | optional (`seq`) | `>=0.3` | `uv sync --extra seq` |

Versions were read with `micromamba run -n rogen-tools <tool> --version`. `tabix` and `bgzip` come from the `htslib=1.24` package. `modkit` is the Bioconda package `ont-modkit=0.6.4`. One create command installed the set together:

```bash
micromamba create -y -n rogen-tools -c conda-forge -c bioconda \
  samtools=1.24 bcftools=1.24 htslib=1.24 minimap2=2.31 ont-modkit=0.6.4
```

## October 2026

Checked on 2026-10-28 for the HG002 short-read dry run
(`scripts/dryrun/hg002_shortread_dryrun.py`) and the ONT versus bisulfite
comparison (`docs/METHYLATION_CLOCK_VALIDATION.md`; `modkit` 0.6.4 in
`rogen-tools`).

| Tool | Found | Version | Install notes |
|------|-------|---------|---------------|
| rtg-tools | yes (repo `.tools/rtg`) | 3.13 (nojre) | Bioconda install was blocked from this host (`conda.anaconda.org` 403). Unpacked `rtg-tools-3.13-nojre.zip` under `.tools/rtg-tools/` with Temurin 17 under `.tools/jdk-17/`. Wrapper: `.tools/rtg`. Prefer `micromamba install -n rogen-tools -c bioconda rtg-tools` when bioconda is reachable. |
| bwa-mem2 | no | — | Optional; dry run default `ALIGNER=minimap2`. Install with `micromamba install -n rogen-tools -c bioconda bwa-mem2` when needed. |
| openjdk / Temurin | yes (`.tools/jdk-17`) | 17.0.13 | Required for RTG nojre build. |

`rtg.cfg` in `.tools/rtg-tools/` sets `RTG_TALKBACK=false` and `RTG_USAGE=false`.
