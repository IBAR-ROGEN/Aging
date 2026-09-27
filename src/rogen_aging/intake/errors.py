"""Errors for the sequencing intake layer.

This package does not make network calls.
"""

from __future__ import annotations


class IntakeError(Exception):
    """Intake refused to continue.

    Callers must not write the processing output or the linkage table after this
    exception. Messages must not include the pseudonym key.
    """
