#!/usr/bin/env python3
import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn("scripts/alphagenome_sequence_comparer.py moved to scripts/alphagenome/alphagenome_sequence_comparer.py", DeprecationWarning, stacklevel=1)
    runpy.run_path(str(Path(__file__).resolve().parent / "alphagenome" / "alphagenome_sequence_comparer.py"), run_name="__main__")
