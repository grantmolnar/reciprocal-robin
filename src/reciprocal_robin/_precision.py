"""Scoped precision for the mpmath point and interval contexts.

Importing :mod:`reciprocal_robin` leaves both contexts unchanged.
:func:`working_precision` changes ``mp.dps`` and ``iv.dps`` together and
restores both on exit. Since mpmath precision is process-global mutable state,
overlapping precision contexts should not be used concurrently.
"""

from contextlib import contextmanager
from typing import Iterator

from mpmath import iv, mp

from ._validation import require_precision

DEFAULT_DPS = 50


def set_default_precision(dps: int = DEFAULT_DPS) -> None:
    """Set both global mpmath contexts to ``dps`` until changed again.

    Library and proof code should normally use :func:`working_precision`, which
    restores the caller's precision on exit.
    """
    precision = require_precision(dps)
    mp.dps = precision
    iv.dps = precision


@contextmanager
def working_precision(dps: int = DEFAULT_DPS) -> Iterator[None]:
    """Temporarily use ``dps`` decimal digits in both mpmath contexts."""
    precision = require_precision(dps)
    prior_mp, prior_iv = mp.dps, iv.dps
    mp.dps = precision
    iv.dps = precision
    try:
        yield
    finally:
        mp.dps = prior_mp
        iv.dps = prior_iv
