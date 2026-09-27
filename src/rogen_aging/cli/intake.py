"""``rogen-intake`` console entry."""

from __future__ import annotations

from pathlib import Path

import typer
from loguru import logger

from rogen_aging.config.cli import config_option, load_cli_config
from rogen_aging.intake.errors import IntakeError
from rogen_aging.intake.key import make_key
from rogen_aging.intake.run import run_intake

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Checksum, pseudonymize, and leak-check a sequencing delivery.",
)


@app.command("run")
def run_cmd(
    manifest: Path = typer.Option(
        ...,
        "--manifest",
        exists=True,
        dir_okay=False,
        help="Delivery TSV: original_sample_id, file_name, file_type, checksum_type, checksum.",
    ),
    delivery_dir: Path = typer.Option(
        ...,
        "--delivery-dir",
        exists=True,
        file_okay=False,
        help="Directory containing the files named in the manifest.",
    ),
    out_dir: Path = typer.Option(
        ...,
        "--out-dir",
        help="Processing directory for pseudonymized files and intake_report.json.",
    ),
    custodian_dir: Path = typer.Option(
        ...,
        "--custodian-dir",
        help="Directory outside the repo and outside --out-dir for linkage.tsv.",
    ),
    config: Path | None = config_option(),
) -> None:
    """Verify checksums, pseudonymize files, and write the linkage table."""
    load_cli_config(config)
    try:
        destination = run_intake(manifest, delivery_dir, out_dir, custodian_dir)
    except IntakeError as exc:
        logger.error("{message}", message=str(exc))
        raise typer.Exit(code=1) from exc
    typer.echo(f"Wrote intake output: {destination}")


@app.command("make-key")
def make_key_cmd(
    out: Path = typer.Option(
        ...,
        "--out",
        help="Destination for a random 32-byte key (mode 600, outside the repository).",
    ),
) -> None:
    """Write a pseudonym key. The key bytes are not printed."""
    try:
        make_key(out)
    except IntakeError as exc:
        logger.error("{message}", message=str(exc))
        raise typer.Exit(code=1) from exc
    typer.echo(f"Wrote key file: {out}")


def main() -> None:
    """Console entry for ``rogen-intake``."""
    app()


if __name__ == "__main__":
    main()
