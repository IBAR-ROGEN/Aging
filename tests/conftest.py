"""Pytest configuration for headless-safe matplotlib and isolated pipeline config."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

MOCK_CONFIG_PATH = Path(__file__).resolve().parent / "fixtures" / "mock_config.yaml"


@pytest.fixture(autouse=True)
def _isolated_pipeline_config() -> None:
    """Load the mock YAML config for every test to isolate from production defaults.

    Uses ``tests/fixtures/mock_config.yaml`` merged on top of ``config/default.yaml``.
    """
    from rogen_aging.config import load_config, reset_config

    load_config(MOCK_CONFIG_PATH)
    yield
    reset_config()


@pytest.fixture
def mock_config_path() -> Path:
    """Absolute path to the pytest mock configuration YAML."""
    return MOCK_CONFIG_PATH


@pytest.fixture
def pipeline_config(mock_config_path: Path):
    """Return the active DictConfig loaded from the mock YAML."""
    from rogen_aging.config import load_config

    return load_config(mock_config_path)
