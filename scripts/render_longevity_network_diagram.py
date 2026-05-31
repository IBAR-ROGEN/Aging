#!/usr/bin/env python3
import runpy, warnings
from pathlib import Path
if __name__ == "__main__":
    warnings.warn("scripts/render_longevity_network_diagram.py moved to scripts/figures/render_longevity_network_diagram.py", DeprecationWarning, stacklevel=1)
    runpy.run_path(str(Path(__file__).resolve().parent / "figures" / "render_longevity_network_diagram.py"), run_name="__main__")
