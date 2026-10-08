"""Reciprocal kappa-CA numbers and their threshold parameters.

``reciprocal_ca(k, eps)`` implements the reciprocal-CA structure theorem: the
local thresholds ``F^[k](p, ell)`` determine the exponent of each prime. The
result is returned as a :class:`Factorization`.
"""

from mpmath import ceil, mpf
from sympy import primerange

from .arithmetic import F_local
from .factorization import Factorization
from ._field import current_field
from ._types import NumericInput, Scalar
from ._validation import (
    require_epsilon_in_ca_domain,
    require_integer,
    require_kappa,
    require_real_greater_than,
)


_MAX_THRESHOLD = mpf(10) ** 20
_BISECTION_ITERATIONS = 120


class ThresholdTooLarge(OverflowError):
    """Raised when a threshold root exceeds the configured search ceiling."""

    def __init__(
        self, k: NumericInput, ell: int, eps: NumericInput, max_threshold: NumericInput
    ) -> None:
        super().__init__(
            "x_threshold: threshold exceeds "
            f"{max_threshold} at k={k}, ell={ell}, eps={eps}; "
            "move eps farther from -1/k or raise max_threshold explicitly"
        )
        self.k = k
        self.ell = ell
        self.eps = eps
        self.max_threshold = max_threshold


class MaxExponentExceeded(ValueError):
    """Raised when ``reciprocal_ca`` reaches its configured exponent cap."""

    def __init__(self, k: NumericInput, eps: NumericInput, max_exponent: int) -> None:
        super().__init__(
            f"reciprocal_ca: max_exponent={max_exponent} reached at "
            f"k={k}, eps={eps}. Raise max_exponent or move eps farther from -1/k."
        )
        self.k = k
        self.eps = eps
        self.max_exponent = max_exponent


def x_threshold(
    k: NumericInput,
    ell: int,
    eps: NumericInput,
    *,
    max_threshold: NumericInput = _MAX_THRESHOLD,
    bisection_iterations: int = _BISECTION_ITERATIONS,
) -> Scalar | None:
    """Solve ``F^[k](x, ell) = eps`` for ``x > 1`` by point bisection.

    Return ``None`` when the root lies below 2. Raise
    :class:`ThresholdTooLarge` when it lies above ``max_threshold``. The search
    uses point arithmetic because interval values cannot order the bisection
    branches.
    """
    require_kappa(k)
    require_epsilon_in_ca_domain(k, eps)
    require_integer(ell, "ell", minimum=1)
    require_real_greater_than(max_threshold, 2, "max_threshold")
    require_integer(bisection_iterations, "bisection_iterations", minimum=1)
    F = current_field()
    if F.name != "mpf":
        raise TypeError("x_threshold is a point bisection; use MpfField explicitly")

    k = F.number(k)
    eps = F.number(eps)
    max_threshold = F.number(max_threshold)
    if F_local(2, ell, k) < eps:
        return None
    lo, hi = F.number(2), F.number(4)
    while F_local(hi, ell, k) > eps:
        hi *= 2
        if hi > max_threshold:
            raise ThresholdTooLarge(k, ell, eps, max_threshold)
    for _ in range(bisection_iterations):
        mid = (lo + hi) / 2
        if F_local(mid, ell, k) > eps:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _to_int_upper(x: Scalar) -> int:
    """Return the least integer greater than or equal to a point value."""
    return int(ceil(x))


def reciprocal_ca(k: NumericInput, eps: NumericInput, max_exponent: int = 200) -> Factorization:
    r"""Return the reciprocal ``k``-CA number ``N^[k](eps)``.

    The reciprocal-CA structure theorem determines the factorization from one
    real parameter ``eps``. At exponent ``ell``, the relevant primes lie in
    ``(x_{ell+1}, x_ell]``, where ``F^[k](x_ell, ell) = eps``.

    The first threshold bounds the candidate primes. For each such prime, the
    largest exponent satisfying ``F^[k](p, ell) >= eps`` is chosen directly;
    equality therefore selects the larger CA integer.

    Parameters
    ----------
    k : numeric
        Real exponent ``> 1``.
    eps : numeric
        Threshold parameter in ``(-1/k, infty)``.
    max_exponent : int
        Finite cap on the exponent search. Raises :class:`MaxExponentExceeded`
        if exceeded; the message reports the active cap. Thresholds
        beyond ``10^20`` raise :class:`ThresholdTooLarge` rather than being
        silently replaced by a finite sentinel.

    Returns
    -------
    Factorization
        The CA number, ready for ``.omega``, ``.p_max``, ``.log_n()``,
        ``.rho_k(k)``, ``.robin_ratio(k)``.

    Examples
    --------
    >>> from mpmath import mpf
    >>> from reciprocal_robin import MpfField, using_field
    >>> from reciprocal_robin.exploration import reciprocal_ca
    >>> with using_field(MpfField()):
    ...     fac = reciprocal_ca(mpf(2), mpf("-0.40"))
    >>> fac.omega, fac.p_max
    (11, 31)
    """
    require_kappa(k)
    require_epsilon_in_ca_domain(k, eps)
    require_integer(max_exponent, "max_exponent", minimum=1)
    F = current_field()
    if F.name != "mpf":
        raise TypeError("reciprocal_ca is a point constructor; use MpfField explicitly")
    k = F.number(k)
    eps = F.number(eps)

    # The first threshold only bounds the candidate primes. Exponents are chosen
    # from the local inequalities themselves, so equality selects the larger
    # local maximizer rather than depending on an approximate root.
    x_one = x_threshold(k, 1, eps)
    if x_one is None or x_one < 2:
        return Factorization._from_trusted_dict({})

    factorization: dict[int, int] = {}
    prime_limit = _to_int_upper(x_one)
    for p in primerange(2, prime_limit + 1):
        exponent = 0
        for ell in range(1, max_exponent + 1):
            if F_local(p, ell, k) >= eps:
                exponent = ell
            else:
                break
        else:
            raise MaxExponentExceeded(k, eps, max_exponent)
        if exponent:
            factorization[p] = exponent
    return Factorization._from_trusted_dict(factorization)
