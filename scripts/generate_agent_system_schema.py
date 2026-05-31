#!/usr/bin/env python3
import runpy, warnings
from pathlib import Path
if __name__ == "__main__":
    warnings.warn("scripts/generate_agent_system_schema.py moved to scripts/figures/generate_agent_system_schema.py", DeprecationWarning, stacklevel=1)
    runpy.run_path(str(Path(__file__).resolve().parent / "figures" / "generate_agent_system_schema.py"), run_name="__main__")
