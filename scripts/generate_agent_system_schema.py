"""Generate Figure 4: Agent System Schema architecture diagram.

This script creates a visualization of the agent system architecture showing:
- User interaction flow
- Cursor IDE and Model Context Protocol integration
- LongevityForest cluster with BioMART, AlphaFold, and STRING databases
"""

from diagrams import Diagram, Cluster, Edge, Node
from diagrams.onprem.client import User
from pathlib import Path
from typing import Optional
import os
import subprocess
import shutil


def find_graphviz() -> Optional[str]:
    """Try to find Graphviz dot executable in common locations."""
    # Check if dot is already in PATH
    dot_path = shutil.which("dot")
    if dot_path:
        return dot_path
    
    # Common installation paths
    common_paths = [
        "/usr/local/bin/dot",
        "/opt/homebrew/bin/dot",
        "/usr/bin/dot",
        "/opt/local/bin/dot",  # MacPorts
        "/Applications/Graphviz.app/Contents/MacOS/dot",
    ]
    
    for path in common_paths:
        if os.path.exists(path) and os.access(path, os.X_OK):
            return path
    
    return None


def create_agent_system_schema(output_path: Optional[str] = None) -> None:
    """Create Figure 4: Agent System Schema architecture diagram.
    
    Args:
        output_path: Path to save the diagram. If None, saves to analysis/ directory.
    """
    # Check for Graphviz
    dot_path = find_graphviz()
    if not dot_path:
        # Try to add common paths to PATH
        common_bin_dirs = [
            "/usr/local/bin",
            "/opt/homebrew/bin",
            "/opt/local/bin",
            "/Applications/Graphviz.app/Contents/MacOS",
        ]
        current_path = os.environ.get("PATH", "")
        for bin_dir in common_bin_dirs:
            if os.path.exists(bin_dir) and bin_dir not in current_path:
                os.environ["PATH"] = f"{bin_dir}:{current_path}"
        
        # Check again
        dot_path = find_graphviz()
        if not dot_path:
            print("Warning: Graphviz not found. Using matplotlib fallback to generate PNG.")
            print("To use Graphviz (better quality), install with: brew install graphviz")
            # Use matplotlib fallback
            import sys
            from pathlib import Path
            fallback_script = Path(__file__).parent / "generate_agent_system_schema_fallback.py"
            if fallback_script.exists():
                import importlib.util
                spec = importlib.util.spec_from_file_location("fallback_module", fallback_script)
                fallback_module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(fallback_module)
                fallback_module.create_agent_system_schema_matplotlib(output_path)
                return
            else:
                raise RuntimeError(
                    "Graphviz not found and fallback script not available.\n"
                    "Please install Graphviz: brew install graphviz\n"
                    "Or download from: https://graphviz.org/download/"
                )
    
    if output_path is None:
        output_dir = Path(__file__).parent.parent / "analysis"
        output_dir.mkdir(exist_ok=True)
        output_path = output_dir / "Fig4_Agent_System_Schema.png"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Convert to absolute path
    output_path_abs = output_path.resolve()
    output_dir_abs = output_path_abs.parent
    output_filename = output_path_abs.stem  # filename without extension
    
    # Change to output directory so diagrams saves there
    original_cwd = os.getcwd()
    try:
        os.chdir(str(output_dir_abs))
        
        # Set graph attributes for a professional look
        graph_attr = {
            "fontsize": "20",
            "bgcolor": "white",
            "pad": "0.5"
        }
        
        with Diagram(
            "Fig 4: LongevityForest Multi-Agent Architecture",
            filename=output_filename,
            show=False,
            direction="LR",
            graph_attr=graph_attr
        ):
            # 1. The User
            user = User("Researcher")

            # 2. The Development Environment
            # We use a generic node for Cursor IDE (VSCode fork)
            ide = Node("Cursor IDE\n(User Interface)", **{"fillcolor": "#007ACC", "style": "rounded,filled"})

            # 3. The Bridge
            # We use a generic node to represent the Model Context Protocol
            mcp = Node("Model Context\nProtocol (MCP)", **{"fillcolor": "#4ECDC4", "style": "rounded,filled"})

            # 4. The Agent Cluster
            with Cluster("LongevityForest Agent Swarm"):
                # We use distinct nodes to represent the type of data/work
                biomart = Node("BioMART Agent\n(Annotations)", **{"fillcolor": "#FFE66D", "style": "rounded,filled"})
                alpha = Node("AlphaFold Agent\n(Structure)", **{"fillcolor": "#95E1D3", "style": "rounded,filled"})
                string = Node("STRING Agent\n(Interactions)", **{"fillcolor": "#FF6B6B", "style": "rounded,filled"})
                
                agents = [biomart, alpha, string]

            # 5. Define the Connections
            # The '<< >>' operator creates a bi-directional edge
            user >> Edge(label="Input", color="black") >> ide
            ide << Edge(label="JSON-RPC", style="bold", color="darkblue") >> mcp
            
            # Connect MCP to all agents
            for agent in agents:
                mcp << Edge(color="darkgreen") >> agent
            
    finally:
        os.chdir(original_cwd)
    
    print(f"Agent System Schema diagram saved to: {output_path_abs}")


if __name__ == "__main__":
    create_agent_system_schema()
