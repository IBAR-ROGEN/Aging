#!/usr/bin/env python3
import runpy, warnings
from pathlib import Path
if __name__ == "__main__":
    warnings.warn("scripts/render_figure1c_mechanisms_network.py moved to scripts/figures/render_figure1c_mechanisms_network.py", DeprecationWarning, stacklevel=1)
    runpy.run_path(str(Path(__file__).resolve().parent / "figures" / "render_figure1c_mechanisms_network.py"), run_name="__main__")
