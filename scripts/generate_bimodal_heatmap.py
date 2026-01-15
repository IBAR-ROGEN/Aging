#!/usr/bin/env python3
"""Generate the bimodal risk heatmap (Figure 2).

This script generates the heatmap described in Activity 2.1.7, showing the 
protective vs. risk effects of candidate longevity genes across different conditions.

Usage:
    python scripts/generate_bimodal_heatmap.py
    or
    uv run python scripts/generate_bimodal_heatmap.py
"""

from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.rogen_aging.methylation_visualizations import create_bimodal_risk_heatmap

if __name__ == "__main__":
    create_bimodal_risk_heatmap()
