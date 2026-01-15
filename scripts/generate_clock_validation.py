#!/usr/bin/env python3
"""Generate the methylation clock validation plot (Figure 3).

This script generates the scatter plot for Activity 2.1.10, showing the
relationship between chronological age and DNAm predicted age with MAE ~ 2.1 years.

Usage:
    python scripts/generate_clock_validation.py
    or
    uv run python scripts/generate_clock_validation.py
"""

from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.rogen_aging.methylation_visualizations import create_clock_validation_plot

if __name__ == "__main__":
    create_clock_validation_plot()
