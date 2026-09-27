# Toolchain inventory — September 2026

Checked on 2026-09-27 on this machine: macOS 27.0, Apple Silicon (`arm64`).

CLI tools were probed with `command -v` and, when present, `--version`. Nothing was installed. Homebrew is not on `PATH` (`/opt/homebrew/bin/brew` and `/usr/local/bin/brew` are absent), and conda, mamba, micromamba, and pixi are absent too. The install column is the command for this OS, not a command that was run.

Python packages were checked with `importlib.metadata` in the project environment (`uv run python`, `.venv`, Python 3.12.8) and in system `python3`. Both environments lack `pysam` and `pod5`.

| Tool | Found | Version | Install command (this OS) |
|------|-------|---------|---------------------------|
| samtools | no | — | `brew install samtools` |
| bcftools | no | — | `brew install bcftools` |
| minimap2 | no | — | `brew install minimap2` |
| modkit | no | — | `conda install bioconda::ont-modkit` |
| tabix | no | — | `brew install htslib` |
| bgzip | no | — | `brew install htslib` |
| pysam | no | — | `uv add pysam` |
| pod5 | no | — | `uv add pod5` |

`tabix` and `bgzip` are HTSlib binaries; Homebrew ships them in the `htslib` formula, not as separate formulae. `samtools` and `bcftools` depend on that formula.

`modkit` has no Homebrew formula. Bioconda package `ont-modkit` publishes an `osx-arm64` build. Upstream GitHub release binaries are Linux-only; a Mac build is from source (`cargo` or the Apple Silicon compile script in [nanoporetech/modkit](https://github.com/nanoporetech/modkit)).
