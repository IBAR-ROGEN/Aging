#!/usr/bin/env python3
"""Pre-flight validate a DNA methylation beta matrix for aging-clock use.

Thin wrapper around ``rogen-clock validate-matrix``. Prefer the console entry:

    uv run rogen-clock validate-matrix \\
      --matrix data/methylation_matrix.csv \\
      --metadata data/metadata.csv \\
      --expected-cpgs data/expected_cpgs.txt \\
      --missing-policy report \\
      --log-out outputs/methylation_preflight.log
"""

from __future__ import annotations

from rogen_aging.cli.clock import app

if __name__ == "__main__":
    # Invoke only the validate-matrix subcommand when this script is run directly
    # with its historical argv shape; otherwise fall through to the full app.
    import sys

    if len(sys.argv) == 1 or sys.argv[1].startswith("-"):
        sys.argv.insert(1, "validate-matrix")
    app()
