# Installing Graphviz for Diagram Generation

To generate PNG files from the diagram scripts, you need to install Graphviz.

## macOS Installation

### Option 1: Homebrew (Recommended)
```bash
brew install graphviz
```

### Option 2: MacPorts
```bash
sudo port install graphviz
```

### Option 3: Conda/Mamba
```bash
conda install -c conda-forge graphviz
# or
mamba install -c conda-forge graphviz
```

### Option 4: Direct Download
Download the installer from: https://graphviz.org/download/

## Verify Installation

After installing, verify Graphviz is available:
```bash
dot -V
```

You should see version information if installed correctly.

## Generate the Diagram

Once Graphviz is installed, run:
```bash
uv run python scripts/generate_agent_system_schema.py
```

The PNG will be saved to `analysis/Fig4_Agent_System_Schema.png`.
