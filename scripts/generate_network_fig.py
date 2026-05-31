#!/usr/bin/env python3
import runpy, warnings
from pathlib import Path
if __name__ == "__main__":
    warnings.warn("scripts/generate_network_fig.py moved to scripts/figures/generate_network_fig.py", DeprecationWarning, stacklevel=1)
    runpy.run_path(str(Path(__file__).resolve().parent / "figures" / "generate_network_fig.py"), run_name="__main__")
