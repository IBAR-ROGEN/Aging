#!/usr/bin/env python3
import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn("scripts/visualize_alphagenome_results.py moved to scripts/alphagenome/visualize_alphagenome_results.py", DeprecationWarning, stacklevel=1)
    runpy.run_path(str(Path(__file__).resolve().parent / "alphagenome" / "visualize_alphagenome_results.py"), run_name="__main__")
