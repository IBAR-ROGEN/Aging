"""``rogen-ukb-integrate`` console entry."""

from __future__ import annotations

from rogen_aging.integration.run_cli import main


def entry() -> None:
    raise SystemExit(main())
