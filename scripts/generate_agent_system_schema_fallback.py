"""Fallback version that uses matplotlib if Graphviz is not available."""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from pathlib import Path
from typing import Optional


def create_agent_system_schema_matplotlib(output_path: Optional[str] = None) -> None:
    """Create Figure 4 using matplotlib as fallback."""
    if output_path is None:
        output_dir = Path(__file__).parent.parent / "analysis"
        output_dir.mkdir(exist_ok=True)
        output_path = output_dir / "Fig4_Agent_System_Schema.png"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(16, 10))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis('off')
    ax.set_title("Fig 4: LongevityForest Multi-Agent Architecture", fontsize=20, pad=20, fontweight='bold')
    
    # Colors
    colors = {
        'user': '#E8F4F8',
        'ide': '#007ACC',
        'mcp': '#4ECDC4',
        'biomart': '#FFE66D',
        'alpha': '#95E1D3',
        'string': '#FF6B6B',
        'cluster': '#F0F0F0',
    }
    
    # 1. User (left)
    user_box = FancyBboxPatch(
        (0.5, 3.5), 1.5, 1,
        boxstyle="round,pad=0.1",
        facecolor=colors['user'],
        edgecolor='black',
        linewidth=2
    )
    ax.add_patch(user_box)
    ax.text(1.25, 4, 'Researcher', ha='center', va='center', fontsize=12, fontweight='bold')
    
    # 2. Cursor IDE (middle-left)
    ide_box = FancyBboxPatch(
        (2.5, 3.5), 1.5, 1,
        boxstyle="round,pad=0.1",
        facecolor=colors['ide'],
        edgecolor='black',
        linewidth=2
    )
    ax.add_patch(ide_box)
    ax.text(3.25, 4, 'Cursor IDE\n(User Interface)', ha='center', va='center', fontsize=11, fontweight='bold', color='white')
    
    # 3. MCP (middle)
    mcp_box = FancyBboxPatch(
        (4.5, 3.5), 1.5, 1,
        boxstyle="round,pad=0.1",
        facecolor=colors['mcp'],
        edgecolor='black',
        linewidth=2
    )
    ax.add_patch(mcp_box)
    ax.text(5.25, 4, 'Model Context\nProtocol (MCP)', ha='center', va='center', fontsize=11, fontweight='bold')
    
    # 4. LongevityForest Cluster (right)
    cluster_box = FancyBboxPatch(
        (6.8, 1.5), 4.5, 5,
        boxstyle="round,pad=0.2",
        facecolor=colors['cluster'],
        edgecolor='black',
        linewidth=2,
        linestyle='--'
    )
    ax.add_patch(cluster_box)
    ax.text(9.05, 6.2, 'LongevityForest Agent Swarm', ha='center', va='center', 
            fontsize=14, fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black'))
    
    # Agents in cluster
    biomart_box = FancyBboxPatch(
        (7.5, 4.5), 1.2, 0.8,
        boxstyle="round,pad=0.1",
        facecolor=colors['biomart'],
        edgecolor='black',
        linewidth=1.5
    )
    ax.add_patch(biomart_box)
    ax.text(8.1, 4.9, 'BioMART Agent\n(Annotations)', ha='center', va='center', fontsize=9, fontweight='bold')
    
    alpha_box = FancyBboxPatch(
        (9.0, 4.5), 1.2, 0.8,
        boxstyle="round,pad=0.1",
        facecolor=colors['alpha'],
        edgecolor='black',
        linewidth=1.5
    )
    ax.add_patch(alpha_box)
    ax.text(9.6, 4.9, 'AlphaFold Agent\n(Structure)', ha='center', va='center', fontsize=9, fontweight='bold')
    
    string_box = FancyBboxPatch(
        (10.5, 4.5), 1.2, 0.8,
        boxstyle="round,pad=0.1",
        facecolor=colors['string'],
        edgecolor='black',
        linewidth=1.5
    )
    ax.add_patch(string_box)
    ax.text(11.1, 4.9, 'STRING Agent\n(Interactions)', ha='center', va='center', fontsize=9, fontweight='bold')
    
    # Connections
    # User -> IDE
    arrow1 = FancyArrowPatch((2.0, 4), (2.5, 4), 
                            arrowstyle='->', lw=2, color='black',
                            mutation_scale=20)
    ax.add_patch(arrow1)
    ax.text(2.25, 4.3, 'Input', ha='center', fontsize=9, fontweight='bold')
    
    # IDE <-> MCP (bidirectional)
    arrow2a = FancyArrowPatch((4.0, 4), (4.5, 4),
                             arrowstyle='->', lw=2.5, color='darkblue',
                             mutation_scale=20)
    ax.add_patch(arrow2a)
    arrow2b = FancyArrowPatch((5.0, 4), (4.5, 4),
                             arrowstyle='->', lw=2, color='darkblue', linestyle='--',
                             mutation_scale=20)
    ax.add_patch(arrow2b)
    ax.text(4.5, 4.3, 'JSON-RPC', ha='center', fontsize=9, fontweight='bold', color='darkblue')
    
    # MCP -> Agents (bidirectional)
    agent_positions = [8.1, 9.6, 11.1]  # x positions of agents
    for i, pos in enumerate(agent_positions):
        # Forward arrow
        arrow_fwd = FancyArrowPatch((6.0, 4), (pos, 4.9),
                                   arrowstyle='->', lw=2, color='darkgreen',
                                   mutation_scale=15)
        ax.add_patch(arrow_fwd)
        # Backward arrow
        arrow_bwd = FancyArrowPatch((pos, 4.9), (6.0, 4),
                                   arrowstyle='->', lw=1.5, color='darkgreen', linestyle='--',
                                   mutation_scale=15)
        ax.add_patch(arrow_bwd)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Agent System Schema diagram saved to: {output_path}")


if __name__ == "__main__":
    create_agent_system_schema_matplotlib()
