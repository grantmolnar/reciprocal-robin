"""Immutable prime factorizations and their normalized arithmetic.

:class:`Factorization` stores canonical prime-exponent pairs and evaluates the
quantities that can be computed directly from an already-factored integer.
"""

from dataclasses import dataclass, field
from typing import Mapping

from ._field import current_field
from ._types import NumericInput, Scalar
from ._validation import require_kappa, validated_factorization
from .arithmetic import rksigma_pp


@dataclass(frozen=True)
class Factorization:
    r"""Immutable prime factorization of a positive integer.

    Instances are hashable, frozen, and ordered by increasing prime. Construct
    them with :meth:`from_dict`; direct construction must already use the same
    canonical form.

    Construction
    ------------
    >>> fac = Factorization.from_dict({2: 2, 3: 2})  # n = 36
    >>> fac.omega
    2
    >>> fac.p_max
    3

    Computation
    -----------
    Numeric methods use :func:`._field.current_field`. The default is interval
    arithmetic; :class:`._field.MpfField` gives point values.

    >>> from mpmath import mpf
    >>> from reciprocal_robin import MpfField, using_field
    >>> with using_field(MpfField()):
    ...     value = fac.rho_k(mpf(2))
    >>> float(value)
    2.2569444444444446
    """

    _pairs: tuple[tuple[int, int], ...] = field()

    def __post_init__(self) -> None:
        """Require direct construction to use canonical prime-exponent pairs."""
        try:
            mapping = dict(self._pairs)
        except (TypeError, ValueError) as exc:
            raise TypeError("_pairs must be an iterable of (prime, exponent) pairs") from exc
        if len(mapping) != len(self._pairs):
            raise ValueError("factorization contains a repeated prime")
        if validated_factorization(mapping) != self._pairs:
            raise ValueError("factorization pairs must be canonical and ordered by prime")

    @classmethod
    def _from_canonical_pairs(cls, pairs: tuple[tuple[int, int], ...]) -> "Factorization":
        """Construct from canonical pairs that have already been validated."""
        instance = object.__new__(cls)
        object.__setattr__(instance, "_pairs", pairs)
        return instance

    @classmethod
    def from_dict(cls, d: Mapping[int, int]) -> "Factorization":
        """Construct from a mapping after validating primes and exponents."""
        return cls._from_canonical_pairs(validated_factorization(d))

    @classmethod
    def _from_trusted_dict(cls, d: Mapping[int, int]) -> "Factorization":
        """Construct from trusted internal prime-exponent data."""
        return cls._from_canonical_pairs(tuple(sorted(d.items())))

    @property
    def primes_to_exponents(self) -> dict[int, int]:
        """Return a mutable copy of the prime-to-exponent mapping."""
        return dict(self._pairs)

    @property
    def omega(self) -> int:
        """Number of distinct prime divisors."""
        return len(self._pairs)

    @property
    def p_max(self) -> int:
        """Largest prime divisor, or 0 for the factorization of 1."""
        return self._pairs[-1][0] if self._pairs else 0

    def log_n(self) -> Scalar:
        """Natural logarithm of the integer represented by this factorization."""
        field = current_field()
        total = field.number(0)
        for p, exponent in self._pairs:
            total = total + exponent * field.log(field.number(p))
        return total

    def rho_k(self, k: NumericInput) -> Scalar:
        """Return ``rho^[k](n)`` for the represented integer ``n``."""
        require_kappa(k)
        field = current_field()
        product = field.number(1)
        for p, exponent in self._pairs:
            product = product * field.number(p) ** exponent * rksigma_pp(p, exponent, k) / k
        return product ** (field.number(1) / (k - 1))

    def robin_ratio(self, k: NumericInput) -> Scalar:
        """Return ``rho^[k](n) / (e^gamma log log n)``."""
        if not self._pairs:
            raise ValueError("robin_ratio is undefined for the factorization of 1")
        require_kappa(k)
        field = current_field()
        log_n = self.log_n()
        if log_n <= field.number(1):
            raise ValueError("robin_ratio requires n >= 3")
        return self.rho_k(k) / (field.exp(field.euler_gamma()) * field.log(log_n))

    def __bool__(self) -> bool:
        return bool(self._pairs)

    def __len__(self) -> int:
        return len(self._pairs)
