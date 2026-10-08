"""The normalized reciprocal-sigma function and its Robin ratio.

Fixed-parameter formulas use the active numerical field. The critical-exponent
solver :func:`kappa_star` is point-valued; boundary classification may instead
be certified with :func:`boundary_difference` under interval arithmetic.
"""

from typing import Mapping, cast

import mpmath
from mpmath import mpf, inf
from sympy import factorint

from .arithmetic import rksigma_pp, sigma_minus_one_pp
from ._field import Trilean, current_field
from ._precision import working_precision
from ._types import NumericInput, Scalar
from ._validation import (
    require_integer,
    require_kappa,
    require_positive,
    require_real_at_least,
)


def rho_k(n: int, k: NumericInput) -> Scalar:
    r"""Return the normalized reciprocal-sigma value ``rho^[k](n)``.

    Defined by

    .. math::
       \rho^{[k]}(n)
       \;=\;
       \left(
         \frac{n \cdot \sigma^{[k]}_{-1}(n)}{k^{\omega(n)}}
       \right)^{1/(k-1)},

    where ``sigma^[k]_{-1} = mu * sigma_{-1}^k`` and ``omega(n)`` counts
    distinct prime divisors. The implementation evaluates this product over the
    prime powers of ``n``.

    Parameters
    ----------
    n : int
        Positive integer.
    k : numeric
        Exponent ``> 1`` accepted by the active numerical field.

    Returns
    -------
    field element
        ``rho^[k](n)`` in the active numerical field.

    Examples
    --------
    >>> from mpmath import mpf
    >>> from reciprocal_robin import MpfField, using_field
    >>> with using_field(MpfField()):
    ...     value = rho_k(36, mpf(2))
    >>> float(value)
    2.2569444444444446
    """
    require_integer(n, "n", minimum=1)
    require_kappa(k)
    F = current_field()
    fac = cast(Mapping[int, int], factorint(n))
    prod = F.number(1)
    for p, a in fac.items():
        prod = prod * F.number(p) ** a * rksigma_pp(p, a, k) / k
    return prod ** (F.number(1) / (k - 1))


def robin_ratio(n: int, k: NumericInput) -> Scalar:
    r"""Normalized Robin ratio ``rho^[k](n) / (e^gamma log log n)``.

    The Robin-type inequality is exactly the assertion that this ratio is
    less than 1.

    Parameters
    ----------
    n : int
        Integer ``>= 3``; ``log log n`` is not defined at ``n=1,2``.
    k : numeric
        Exponent ``> 1`` accepted by the active numerical field.

    Returns
    -------
    field element
        Dimensionless ratio in the active numerical field.
    """
    require_integer(n, "n", minimum=3)
    require_kappa(k)
    F = current_field()
    return rho_k(n, k) / (F.exp(F.euler_gamma()) * F.log(F.log(F.number(n))))


def kappa_star(
    n: int,
    k_lo: NumericInput = mpf("1"),
    k_hi: NumericInput = mpf("1000"),
    tol: NumericInput = mpf("1e-20"),
    k_hi_max: NumericInput = mpf("1e12"),
) -> Scalar:
    r"""Numerical point estimate of the critical exponent ``kappa*(n)``.

    The mathematical trichotomy is handled before bisection:

    - ``kappa*(n) = 1`` when the boundary value ``rho^[1](n)`` is at
      least ``e^gamma log log n``;
    - ``kappa*(n) = +infinity`` when ``sigma_{-1}(n)`` is at most that
      scale;
    - otherwise strict monotonicity gives a unique interior root.

    This is a point-arithmetic root finder, not an interval certificate. It
    therefore requires :class:`MpfField`.

    A returned infinity in the interior-search branch means only that no
    violation was found before ``k_hi_max``.
    """
    F = current_field()
    if F.name != "mpf":
        raise TypeError("kappa_star is a point root finder; use MpfField explicitly")

    require_integer(n, "n", minimum=3)

    one = F.number(1)
    k_lo = F.number(k_lo)
    k_hi = F.number(k_hi)
    tol = F.number(tol)
    k_hi_max = F.number(k_hi_max)
    require_real_at_least(k_lo, 1, "k_lo")
    require_kappa(k_hi, "k_hi")
    require_positive(tol, "tol")
    require_kappa(k_hi_max, "k_hi_max")
    if k_hi_max < k_hi:
        raise ValueError("k_hi_max must be at least k_hi")

    scale = F.exp(F.euler_gamma()) * F.log(F.log(F.number(n)))
    if rho_hat_one(n) >= scale:
        return one

    fac = cast(Mapping[int, int], factorint(n))
    sigma_minus_one = F.number(1)
    for p, a in fac.items():
        sigma_minus_one *= sigma_minus_one_pp(p, a)
    if sigma_minus_one <= scale:
        return inf

    # kappa=1 is a limiting endpoint; every bisection evaluation lies above it.
    k_lo = max(k_lo, one)
    while robin_ratio(n, k_hi) < one:
        if k_hi >= k_hi_max:
            return inf
        k_hi = min(k_hi * 2, k_hi_max)

    while k_hi - k_lo > tol:
        mid = (k_lo + k_hi) / 2
        if robin_ratio(n, mid) < one:
            k_lo = mid
        else:
            k_hi = mid
    return (k_lo + k_hi) / 2


def exceptional_set(k: NumericInput, n_max: int, n_min: int = 3) -> list[int]:
    r"""Enumerate ``E(k) intersected with [n_min, n_max]``.

    Returns the sorted list of integers ``n in [n_min, n_max]`` for
    which ``rho^[k](n) >= e^gamma log log n``, i.e. those that
    *violate* the Robin-type inequality at the given exponent.

    This is a finite brute-force search, with one :func:`robin_ratio`
    evaluation per integer in the requested range.

    Parameters
    ----------
    k : int, float, str, or mpf
        Exponent ``> 1``.
    n_max : int
        Upper bound of the search range, inclusive.
    n_min : int, optional
        Lower bound of the search range, inclusive. Default ``3``
        (since ``robin_ratio`` is undefined for ``n in {1, 2}``).

    Returns
    -------
    list of int
        The sorted exceptional integers in ``[n_min, n_max]``.

    Examples
    --------
    >>> from mpmath import mpf
    >>> exceptional_set(mpf(2), 30)
    [3, 4, 5, 6, 8, 12, 16, 24]
    """
    require_kappa(k)
    require_integer(n_min, "n_min", minimum=3)
    require_integer(n_max, "n_max", minimum=3)
    if n_max < n_min:
        raise ValueError("n_max must be at least n_min")
    k_mpf = mpf(k)
    F = current_field()

    def _ratio_ge_one(nn: int) -> bool:
        r = robin_ratio(nn, k_mpf)
        v = F.is_strictly_less_than_one(r)
        if v == Trilean.AMBIGUOUS:
            with working_precision(mpmath.mp.dps * 4):
                r = robin_ratio(nn, k_mpf)
                v = F.is_strictly_less_than_one(r)
        # An unresolved interval is retained: the finite search must not discard
        # a possible exception.
        return v in (Trilean.NO, Trilean.AMBIGUOUS)

    return [n for n in range(n_min, n_max + 1) if _ratio_ge_one(n)]


def rho_hat_one(n: int) -> Scalar:
    r"""Boundary limit ``rho^[1](n) := lim_{k -> 1+} rho^[k](n)``.

    At a prime power, write

    .. math::
       S_a = \sigma_{-1}(p^a), \qquad u = S_{a-1}/S_a.

    The local power-mean representation gives the exact boundary value

    .. math::
       \rho^{[1]}(p^a)
       = S_a \exp\left(-1-\frac{u\log u}{1-u}\right).

    Multiplying these local factors avoids the less stable equivalent formula
    obtained by differentiating the divisor convolution.
    """
    F = current_field()
    require_integer(n, "n", minimum=1)
    if n == 1:
        return F.number(1)

    result = F.number(1)
    one = F.number(1)
    fac = cast(Mapping[int, int], factorint(n))
    for p, a in fac.items():
        s_a = sigma_minus_one_pp(p, a)
        s_prev = sigma_minus_one_pp(p, a - 1)
        u = s_prev / s_a
        local = s_a * F.exp(-one - u * F.log(u) / (one - u))
        result *= local
    return result


def boundary_difference(n: int) -> Scalar:
    r"""Return ``rho^[1](n) - e^gamma log log n``.

    Under :class:`MpiField`, the returned interval can certify the boundary
    classification whenever it lies strictly on one side of zero.
    """
    F = current_field()
    require_integer(n, "n", minimum=3)
    scale = F.exp(F.euler_gamma()) * F.log(F.log(F.number(n)))
    return rho_hat_one(n) - scale
