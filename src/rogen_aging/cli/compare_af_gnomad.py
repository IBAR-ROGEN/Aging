"""``rogen-compare-af-gnomad`` console entry."""

from __future__ import annotations

from rogen_aging.ukb.gnomad import main


def entry() -> None:
    raise SystemExit(main())


if __name__ == "__main__":
    entry()
