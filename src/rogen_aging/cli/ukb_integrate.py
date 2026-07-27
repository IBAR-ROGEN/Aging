"""``rogen-ukb-integrate`` console entry."""

from __future__ import annotations

from rogen_aging.ukb_integration.run_cli import main


def entry() -> None:
    """Console entry for ``rogen-ukb-integrate``."""
    raise SystemExit(main())
