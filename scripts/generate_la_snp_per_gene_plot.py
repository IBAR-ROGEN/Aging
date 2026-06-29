#!/usr/bin/env python3
import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn("scripts/generate_la_snp_per_gene_plot.py moved to scripts/figures/generate_la_snp_per_gene_plot.py", DeprecationWarning, stacklevel=1)
    runpy.run_path(str(Path(__file__).resolve().parent / "figures" / "generate_la_snp_per_gene_plot.py"), run_name="__main__")
