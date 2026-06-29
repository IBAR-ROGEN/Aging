#!/usr/bin/env python3
import runpy
import warnings
from pathlib import Path

if __name__ == "__main__":
    warnings.warn("scripts/render_dashboard_figure_mockup.py moved to scripts/figures/render_dashboard_figure_mockup.py", DeprecationWarning, stacklevel=1)
    runpy.run_path(str(Path(__file__).resolve().parent / "figures" / "render_dashboard_figure_mockup.py"), run_name="__main__")
