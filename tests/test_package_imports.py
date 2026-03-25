"""Smoke tests: package is importable after installation."""

from __future__ import annotations


def test_import_package() -> None:
    import rogen_aging  # noqa: PLC0415

    assert hasattr(rogen_aging, "__all__")


def test_import_submodules() -> None:
    import rogen_aging.methylation_visualizations  # noqa: PLC0415
    import rogen_aging.network_visualizer  # noqa: PLC0415
    import rogen_aging.pipeline  # noqa: PLC0415

    assert rogen_aging.pipeline.__all__ == []
