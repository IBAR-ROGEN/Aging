"""Tests for ``rogen_aging.config`` loading and CLI wiring."""

from __future__ import annotations

from pathlib import Path

import pytest
from omegaconf import OmegaConf

from rogen_aging.config import (
    cfg_path,
    default_config_path,
    get_config,
    load_config,
    production_config_path,
    risk_weights,
    target_tissues,
)
from rogen_aging.config.cli import load_cli_config
from rogen_aging.integrative.phenotype_integrator import DEFAULT_WEIGHTS, PhenotypeIntegrator


def test_default_config_loads() -> None:
    cfg = load_config()
    assert Path(str(cfg.repo.root)).is_dir()
    assert default_config_path().is_file()
    assert float(cfg.integrative.risk_weights.vep_impact) == pytest.approx(0.25)
    assert cfg_path(cfg, "paths", "models", "clock_elasticnet").name.endswith(".pkl")


def test_production_profile_merges() -> None:
    cfg = load_config(profile="production")
    assert production_config_path().is_file()
    assert str(cfg.apis.ensembl_rest).startswith("https://")


def test_mock_config_isolates_paths(mock_config_path: Path) -> None:
    cfg = load_config(mock_config_path)
    model = cfg_path(cfg, "paths", "models", "clock_elasticnet")
    assert "mock_clock" in model.name
    assert (
        cfg_path(cfg, "paths", "integrative", "output_dir")
        .as_posix()
        .endswith("test_data/integrative/results")
    )


def test_override_weights_sync_module_defaults(tmp_path: Path) -> None:
    override = tmp_path / "weights.yaml"
    override.write_text(
        "\n".join(
            [
                "integrative:",
                "  risk_weights:",
                "    vep_impact: 0.40",
                "    alphagenome: 0.20",
                "    alphamissense: 0.20",
                "    gtex_eqtl: 0.10",
                "    epigenetic: 0.10",
                "  target_tissues:",
                "    - Brain_Cortex",
                "    - Whole_Blood",
                "",
            ]
        ),
        encoding="utf-8",
    )
    cfg = load_config(override)
    weights = risk_weights(cfg)
    assert weights["vep_impact"] == pytest.approx(0.40)
    assert target_tissues(cfg) == ("Brain_Cortex", "Whole_Blood")
    assert DEFAULT_WEIGHTS["vep_impact"] == pytest.approx(0.40)
    integrator = PhenotypeIntegrator()
    assert integrator.weights["vep_impact"] == pytest.approx(0.40)


def test_cli_config_option_loads_override(tmp_path: Path) -> None:
    override = tmp_path / "override.yaml"
    override.write_text("clock:\n  random_state: 99\n", encoding="utf-8")
    cfg = load_cli_config(override)
    assert int(cfg.clock.random_state) == 99
    assert OmegaConf.select(get_config(), "clock.random_state") == 99
