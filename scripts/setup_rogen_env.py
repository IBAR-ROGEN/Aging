#!/usr/bin/env python3
"""
Secure Directory Generator (setup_rogen_env.py)

Creates a standardized bioinformatics project folder structure and ensures
.gitignore blocks data/ and results/ folders plus sensitive file types (.vcf, .bam, .csv)
per the rule: keep data and code strictly separated.

See: https://github.com/IBAR-ROGEN/Aging
"""

from pathlib import Path
from typing import Sequence

import typer

app = typer.Typer(help="Setup standardized bioinformatics project structure with secure .gitignore")

# Standard folder structure for bioinformatics projects
DEFAULT_DIRS: Sequence[str] = (
    "src",
    "data/raw",
    "data/processed",
    "results",
    "docs",
)

# .gitignore block for bioinformatics data/code separation
GITIGNORE_BIOINFORMATICS = """# === Bioinformatics: data/code separation (setup_rogen_env.py) ===
# Block entire data and results folders - never commit large/sensitive data
data/
results/
outputs/

# Block sensitive bioinformatics file types globally (safety net)
*.vcf
*.vcf.gz
*.vcf.bgz
*.bam
*.bam.bai
*.sam
*.cram
*.csv
*.tsv
*.bed
*.bed.gz
*.fastq
*.fastq.gz
*.fq
*.fq.gz
"""

GITIGNORE_SECTION_MARKER = "# === Bioinformatics: data/code separation (setup_rogen_env.py) ==="


def create_directories(root: Path, dirs: Sequence[str]) -> list[Path]:
    """Create directory structure. Returns list of created paths."""
    created: list[Path] = []
    for d in dirs:
        path = root / d
        path.mkdir(parents=True, exist_ok=True)
        if path.exists():
            created.append(path)
    return created


def ensure_gitignore(root: Path) -> bool:
    """
    Ensure .gitignore contains the bioinformatics block.
    Appends the block if the section is not present.
    Returns True if changes were made.
    """
    gitignore_path = root / ".gitignore"
    content = gitignore_path.read_text() if gitignore_path.exists() else ""

    if GITIGNORE_SECTION_MARKER in content:
        return False

    # Append the bioinformatics section
    to_append = "\n" + GITIGNORE_BIOINFORMATICS
    gitignore_path.parent.mkdir(parents=True, exist_ok=True)
    with gitignore_path.open("a") as f:
        f.write(to_append)
    return True


@app.command()
def main(
    root: Path = typer.Option(
        Path.cwd(),
        "--root",
        "-r",
        path_type=Path,
        help="Project root directory (default: current directory)",
    ),
) -> None:
    """Create standardized folder structure and update .gitignore."""
    root = root.resolve()
    if not root.is_dir():
        typer.echo(f"Error: {root} is not a directory", err=True)
        raise typer.Exit(1)

    dirs_created = create_directories(root, DEFAULT_DIRS)
    typer.echo(f"Directories ensured under {root}:")
    for p in sorted(dirs_created, key=lambda x: str(x)):
        typer.echo(f"  {p.relative_to(root)}")

    gitignore_updated = ensure_gitignore(root)
    if gitignore_updated:
        typer.echo("\n.gitignore updated with bioinformatics data/code separation rules.")
    else:
        typer.echo("\n.gitignore already contains the bioinformatics section.")

    typer.echo("\nDone. data/, results/, and *.vcf, *.bam, *.csv (etc.) are blocked from git.")


if __name__ == "__main__":
    app()
