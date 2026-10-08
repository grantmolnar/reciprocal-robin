"""Point and interval arithmetic behind a common field interface.

Arithmetic formulas obtain their numerical field from :func:`current_field`.
The caller selects that field with :func:`using_field`, whose context-local
state is restored on exit. :class:`RealField` contains only the operations that
differ between point and interval arithmetic; ordinary arithmetic remains on
the underlying mpmath values.

Interval comparisons use :class:`Trilean`. An interval wholly below the
threshold gives :attr:`Trilean.YES`, one wholly above it gives
:attr:`Trilean.NO`, and an interval crossing it gives
:attr:`Trilean.AMBIGUOUS`.
"""

from contextlib import contextmanager
from contextvars import ContextVar
from enum import Enum
from typing import Iterator, Protocol, runtime_checkable

import mpmath

from ._types import NumericInput, Scalar


class Trilean(Enum):
    """Result of a comparison whose interval may cross the threshold.

    Point arithmetic returns only :attr:`YES` or :attr:`NO`. Interval
    arithmetic returns :attr:`AMBIGUOUS` when the comparison is unresolved at
    the current precision.
    """

    YES = "yes"
    NO = "no"
    AMBIGUOUS = "ambiguous"


@runtime_checkable
class RealField(Protocol):
    """Constants and transcendental operations that depend on field semantics.

    Field elements themselves support Python's ordinary arithmetic operators.
    """

    name: str

    def number(self, x: NumericInput) -> Scalar:
        """Construct a field element from an int, str, or compatible value."""
        ...

    def log(self, x: Scalar) -> Scalar:
        """Natural logarithm."""
        ...

    def exp(self, x: Scalar) -> Scalar:
        """Exponential."""
        ...

    def euler_gamma(self) -> Scalar:
        """The Euler--Mascheroni constant gamma as a field element."""
        ...

    def is_strictly_less_than_one(self, x: Scalar) -> Trilean:
        """Trilean comparison ``x < 1``."""
        ...


class MpfField:
    """Point arithmetic using mpmath ``mpf`` values.

    Precision comes from ``mp.dps``; :func:`._precision.working_precision`
    changes it temporarily. Comparisons are binary.
    """

    name = "mpf"

    def number(self, x: NumericInput) -> Scalar:
        return mpmath.mpf(x)

    def log(self, x: Scalar) -> Scalar:
        return mpmath.log(x)

    def exp(self, x: Scalar) -> Scalar:
        return mpmath.exp(x)

    def euler_gamma(self) -> Scalar:
        return mpmath.euler

    def is_strictly_less_than_one(self, x: Scalar) -> Trilean:
        return Trilean.YES if x < 1 else Trilean.NO


class MpiField:
    """Directed interval arithmetic using mpmath ``mpi`` values.

    Arithmetic and transcendental operations propagate enclosures at the
    current interval precision. :meth:`is_strictly_less_than_one` returns
    :attr:`Trilean.YES` when the upper endpoint is below 1,
    :attr:`Trilean.NO` when the lower endpoint is at least 1, and
    :attr:`Trilean.AMBIGUOUS` otherwise.

    Examples
    --------
    >>> F = MpiField()
    >>> F.is_strictly_less_than_one(F.number("0.5"))
    <Trilean.YES: 'yes'>
    >>> F.is_strictly_less_than_one(F.number("2"))
    <Trilean.NO: 'no'>
    """

    name = "mpi"

    def number(self, x: NumericInput) -> Scalar:
        return mpmath.mpi(x)

    def log(self, x: Scalar) -> Scalar:
        # mpmath.log uses the point context and rejects nondegenerate mpi
        # values; the interval context preserves the enclosure.
        return mpmath.iv.log(x)

    def exp(self, x: Scalar) -> Scalar:
        return mpmath.iv.exp(x)

    def euler_gamma(self) -> Scalar:
        # Coerce the interval-context constant to mpi so it composes with
        # ordinary mpi arithmetic.
        return mpmath.mpi(mpmath.iv.euler)

    def is_strictly_less_than_one(self, x: Scalar) -> Trilean:
        # An mpi has .a (lower endpoint) and .b (upper endpoint).
        if x.b < 1:
            return Trilean.YES
        if x.a >= 1:
            return Trilean.NO
        return Trilean.AMBIGUOUS


# Interval arithmetic is the default so an unresolved comparison remains
# unresolved rather than being rounded to a binary answer. Point-valued work
# selects MpfField explicitly.
_current_field: ContextVar[RealField] = ContextVar("_current_field", default=MpiField())


def current_field() -> RealField:
    """Return the field currently in scope. Defaults to :class:`MpiField`."""
    return _current_field.get()


def _to_float(x: Scalar) -> float:
    """Return a point value for human-readable exploratory output.

    Intervals are replaced by their midpoint because ``float()`` rejects
    nondegenerate ``mpi`` values. Proof decisions must use interval endpoints
    or :class:`Trilean` comparisons instead.
    """
    if hasattr(x, "a") and hasattr(x, "b"):
        return float((mpmath.mpf(x.a) + mpmath.mpf(x.b)) / 2)
    return float(x)


@contextmanager
def using_field(field: RealField) -> Iterator[RealField]:
    """Run a block with ``field`` as the current numerical backend.

    Examples
    --------
    >>> from reciprocal_robin._field import using_field, MpfField, current_field
    >>> with using_field(MpfField()):
    ...     assert current_field().name == "mpf"
    """
    token = _current_field.set(field)
    try:
        yield field
    finally:
        _current_field.reset(token)
