#!/usr/bin/env python3
"""Generate visualizations for the ROGEN methylation pipeline.

This script creates three visualization files:
1. Pipeline workflow diagram
2. Example DMR analysis visualizations
3. Pipeline summary diagram

Usage:
    python scripts/generate_methylation_visualizations.py
    or
    uv run python scripts/generate_methylation_visualizations.py
"""

from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.rogen_aging.methylation_visualizations import generate_all_visualizations

if __name__ == "__main__":
    generate_all_visualizations()
