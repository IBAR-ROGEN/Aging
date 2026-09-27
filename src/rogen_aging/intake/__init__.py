"""Privacy-first sequencing intake.

No Romanian sequencing data is read here. Callers supply a delivery
manifest and files. The pseudonym key stays outside the repository.
"""

from __future__ import annotations

from rogen_aging.intake.errors import IntakeError
from rogen_aging.intake.key import make_key
from rogen_aging.intake.pseudonym import pseudonym_for
from rogen_aging.intake.run import run_intake

__all__ = [
    "IntakeError",
    "make_key",
    "pseudonym_for",
    "run_intake",
]
