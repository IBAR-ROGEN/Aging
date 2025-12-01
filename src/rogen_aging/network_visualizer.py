"""Network Hub Visualizer for Protein Interaction Topology."""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def create_network_visualization(output_path: str = "Network_Analysis_Nov.png") -> None:
    """
    Create and save a network visualization of protein interactions.
    
    Args:
        output_path: Path where the visualization will be saved
    """
    # Genes identified in previous reports (Longevity + Neurodegeneration overlap)
    # We categorize them to color-code the nodes
    genes = {
        'Hubs': ['APOE', 'HSPA1A', 'TP53', 'MAPT', 'AKT1'],
        'Longevity': ['FOXO3', 'SIRT1', 'GPX1', 'SOD2', 'CAT'],
        'Neuro': ['SNCA', 'PINK1', 'LRRK2', 'PARK7', 'BACE1']
    }

    G = nx.Graph()

    # Add nodes with attributes
    for category, gene_list in genes.items():
        for gene in gene_list:
            G.add_node(gene, category=category)

    # Define edges (simulating interactions found via MCP/BioThings setup)
    edges = [
        ('APOE', 'MAPT'), ('APOE', 'BACE1'), ('TP53', 'SIRT1'), ('TP53', 'FOXO3'),
        ('HSPA1A', 'MAPT'), ('HSPA1A', 'SNCA'), ('HSPA1A', 'LRRK2'),
        ('AKT1', 'FOXO3'), ('AKT1', 'TP53'), ('SIRT1', 'FOXO3'),
        ('GPX1', 'SOD2'), ('CAT', 'SOD2'), ('PINK1', 'PARK7')
    ]
    G.add_edges_from(edges)

    # Compute centrality to determine node size
    centrality = nx.degree_centrality(G)

    # Plotting
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G, k=0.5, seed=42)

    # Draw nodes by category
    colors = {'Hubs': '#FF6B6B', 'Longevity': '#4ECDC4', 'Neuro': '#FFE66D'}
    for cat, color in colors.items():
        nodelist = [n for n, d in G.nodes(data=True) if d['category'] == cat]
        nx.draw_networkx_nodes(
            G, pos, nodelist=nodelist, node_color=color,
            node_size=[centrality[n]*3000 + 500 for n in nodelist],
            label=cat, alpha=0.9
        )

    nx.draw_networkx_edges(G, pos, width=2, alpha=0.4, edge_color='gray')
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold')

    plt.title("Protein Interaction Topology: The 'Resilience Core'")
    plt.legend(scatterpoints=1)
    plt.axis('off')
    
    # Save the figure
    output_file = Path(output_path)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Network visualization saved to: {output_file.absolute()}")
    
    plt.show()


if __name__ == "__main__":
    create_network_visualization()

