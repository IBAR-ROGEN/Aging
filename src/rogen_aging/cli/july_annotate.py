"""``rogen-july-annotate`` console entry."""

from __future__ import annotations

from rogen_aging.annotation.july import app


def entry() -> None:
    """Console entry for ``rogen-july-annotate``."""
    app()


if __name__ == "__main__":
    entry()
