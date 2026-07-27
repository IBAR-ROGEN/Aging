"""``rogen-vcf-synthetic`` console entry."""

from __future__ import annotations

from rogen_aging.vcf.synthetic import main


def entry() -> None:
    main()


if __name__ == "__main__":
    entry()
