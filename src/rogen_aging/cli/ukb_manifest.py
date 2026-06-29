"""``rogen-ukb-manifest`` console entry."""

from __future__ import annotations

from rogen_aging.ukb.manifest import main


def entry() -> None:
    raise SystemExit(main())


if __name__ == "__main__":
    entry()
