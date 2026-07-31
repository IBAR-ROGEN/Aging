"""Shared Typer helpers for ``--config`` wiring."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import typer
from omegaconf import DictConfig

from rogen_aging.config.loader import load_config


def config_option() -> Any:
    """Return a Typer option for ``--config`` / ``-c``."""
    return typer.Option(
        None,
        "--config",
        "-c",
        help="YAML config file merged on top of config/default.yaml.",
        exists=False,
        dir_okay=False,
        resolve_path=True,
    )


def load_cli_config(config: Path | None) -> DictConfig:
    """Load the active pipeline config for a CLI invocation.

    Args:
        config: Optional YAML override path from ``--config``.

    Returns:
        Merged OmegaConf ``DictConfig`` (also set as the process-wide active config).
    """
    return load_config(config)
