#!/usr/bin/env python3
import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn("scripts/generate_pipeline_diagram.py moved to scripts/figures/generate_pipeline_diagram.py", DeprecationWarning, stacklevel=1)
    runpy.run_path(str(Path(__file__).resolve().parent / "figures" / "generate_pipeline_diagram.py"), run_name="__main__")
