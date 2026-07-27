"""Smoke tests: package is importable after installation."""

from __future__ import annotations

import warnings


def test_import_package() -> None:
    import rogen_aging  # noqa: PLC0415

    assert hasattr(rogen_aging, "__all__")


def test_import_submodules() -> None:
    import rogen_aging.clock  # noqa: PLC0415
    import rogen_aging.integrative  # noqa: PLC0415
    import rogen_aging.methylation_visualizations  # noqa: PLC0415
    import rogen_aging.network_visualizer  # noqa: PLC0415
    import rogen_aging.ukb  # noqa: PLC0415
    import rogen_aging.ukb_integration  # noqa: PLC0415
    import rogen_aging.vcf  # noqa: PLC0415

    assert "train_clock" in rogen_aging.clock.__all__
    assert "generate_ukb_rap_mock" in rogen_aging.ukb.__all__
    assert "run_integrative_pipeline" in rogen_aging.integrative.__all__
    assert "run_integration_pipeline" in rogen_aging.ukb_integration.__all__


def test_integration_compat_shim_emits_deprecation() -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        import rogen_aging.integration  # noqa: PLC0415

        assert "run_integration_pipeline" in rogen_aging.integration.__all__
    assert any(issubclass(w.category, DeprecationWarning) for w in caught)
