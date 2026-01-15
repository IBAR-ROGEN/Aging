"""Generate a professional pipeline diagram using the diagrams library.

This script creates a visualization of a bioinformatics pipeline showing:
- Data flow from Nanopore pod5 through Dorado and Modkit
- Apache Parquet storage
- Polars and DuckDB querying
- All orchestrated by Dagster
"""

from diagrams import Diagram, Cluster, Edge
from diagrams.generic.storage import Storage
from diagrams.programming.language import Python
from diagrams.onprem.database import Duckdb
from diagrams.generic.database import SQL
from pathlib import Path
from typing import Optional
import os


def create_pipeline_diagram(output_path: Optional[str] = None) -> None:
    """Create a professional pipeline diagram for scientific reports.
    
    Args:
        output_path: Path to save the diagram. If None, saves to analysis/ directory.
    """
    if output_path is None:
        output_dir = Path(__file__).parent.parent / "analysis"
        output_dir.mkdir(exist_ok=True)
        output_path = output_dir / "Bioinformatics_Pipeline_Diagram.png"
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
        
        with Diagram(
            "Bioinformatics Pipeline Architecture",
            filename=output_filename,
            show=False,
            direction="LR",
            graph_attr={
                "bgcolor": "white",
                "pad": "1.5",
                "splines": "ortho",
                "nodesep": "1.0",
                "ranksep": "1.5",
                "fontsize": "16",
                "fontname": "Arial",
                "labeljust": "l",
            },
            node_attr={
                "fontsize": "11",
                "fontname": "Arial",
                "shape": "box",
                "style": "rounded,filled",
                "fillcolor": "lightblue",
            },
        ):
            with Cluster("Dagster Orchestration", graph_attr={
                "bgcolor": "lightgray",
                "style": "rounded,filled",
                "labeljust": "l",
                "fontsize": "14",
                "fontname": "Arial Bold",
            }):
                # Input: Nanopore pod5
                pod5 = Storage(
                    "Nanopore\npod5",
                    **{"fillcolor": "#E8F4F8", "style": "rounded,filled"}
                )
                
                # Processing tools
                dorado = Python(
                    "Dorado\nBasecalling",
                    **{"fillcolor": "#4ECDC4", "style": "rounded,filled"}
                )
                modkit = Python(
                    "Modkit\nMethylation\nCalling",
                    **{"fillcolor": "#4ECDC4", "style": "rounded,filled"}
                )
                
                # Storage: Apache Parquet
                parquet = Storage(
                    "Apache\nParquet\nStorage",
                    **{"fillcolor": "#FFE66D", "style": "rounded,filled"}
                )
                
                # Query engines
                polars = SQL(
                    "Polars\nQuery Engine",
                    **{"fillcolor": "#95E1D3", "style": "rounded,filled"}
                )
                duckdb = Duckdb(
                    "DuckDB\nQuery Engine",
                    **{"fillcolor": "#95E1D3", "style": "rounded,filled"}
                )
                
                # Flow connections - processing pipeline
                pod5 >> Edge(
                    label="Raw reads",
                    style="bold",
                    color="#2C3E50",
                    penwidth="2.5"
                ) >> dorado
                
                dorado >> Edge(
                    label="BAM files",
                    style="bold",
                    color="#27AE60",
                    penwidth="2.5"
                ) >> modkit
                
                modkit >> Edge(
                    label="Methylation calls",
                    style="bold",
                    color="#E74C3C",
                    penwidth="2.5"
                ) >> parquet
                
                # Query connections - parallel querying
                parquet >> Edge(
                    label="Query",
                    style="dashed",
                    color="#7F8C8D",
                    penwidth="2.0"
                ) >> polars
                
                parquet >> Edge(
                    label="Query",
                    style="dashed",
                    color="#7F8C8D",
                    penwidth="2.0"
                ) >> duckdb
    finally:
        os.chdir(original_cwd)
    
    print(f"Pipeline diagram saved to: {output_path_abs}")


if __name__ == "__main__":
    create_pipeline_diagram()
