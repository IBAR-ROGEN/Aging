"""``rogen-ukb-mock-rap`` console entry."""

from __future__ import annotations

from rogen_aging.ukb.mock_rap import app


def entry() -> None:
    app()


if __name__ == "__main__":
    entry()
