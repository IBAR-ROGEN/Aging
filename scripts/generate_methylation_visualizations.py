#!/usr/bin/env python3
import runpy, warnings
from pathlib import Path
if __name__ == "__main__":
    warnings.warn("scripts/generate_methylation_visualizations.py moved to scripts/figures/generate_methylation_visualizations.py", DeprecationWarning, stacklevel=1)
    runpy.run_path(str(Path(__file__).resolve().parent / "figures" / "generate_methylation_visualizations.py"), run_name="__main__")
