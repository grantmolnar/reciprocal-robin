"""Runtime support for certification commands.

Some canonical endpoints have more than 4300 decimal digits, beyond CPython's
default limit for decimal integer conversion. Certification temporarily removes
that limit while serializing or reconstructing those witnesses, then restores
the previous setting. The setting is process-global rather than thread-local.
"""

from __future__ import annotations

import sys
from collections.abc import Iterator
from contextlib import contextmanager


@contextmanager
def unbounded_decimal_integers() -> Iterator[None]:
    """Allow arbitrary-length decimal integer conversion within this context.

    Each invocation restores the setting it found on entry, so nesting is safe.
    Concurrent uses are not isolated because CPython's setting is process-global.
    On interpreters without this guard, the context manager is a no-op.
    """

    getter = getattr(sys, "get_int_max_str_digits", None)
    setter = getattr(sys, "set_int_max_str_digits", None)
    if getter is None or setter is None:
        yield
        return

    previous = getter()
    setter(0)
    try:
        yield
    finally:
        setter(previous)
