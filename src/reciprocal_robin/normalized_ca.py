"""Prime-local formulas for extremizers of ``rho^[k](n) n^{-eps}``.

These are not the reciprocal-CA objects of :mod:`reciprocal_robin.ca`, which
maximize ``sigma_{-1}^{[k]}(n) n^{-k eps}``. The normalization here is the one
used in the Robin endpoint-convexity argument.
"""

from __future__ import annotations

from ._field import current_field
from ._types import NumericInput, Scalar
from .arithmetic import rksigma_pp
from ._validation import require_finite, require_integer, require_kappa, require_prime


def rho_pp(p: int, a: int, k: NumericInput) -> Scalar:
    r"""Return the prime-power factor ``rho^[k](p^a)``.

    By convention ``rho_pp(p, 0, k) = 1``.
    """
    require_prime(p)
    require_integer(a, "a", minimum=0)
    require_kappa(k)
    F = current_field()
    if a == 0:
        return F.number(1)
    p_f = F.number(p)
    k_f = F.number(k)
    base = p_f**a * rksigma_pp(p, a, k_f) / k_f
    return base ** (F.number(1) / (k_f - 1))


def normalized_local_objective(p: int, a: int, k: NumericInput, eps: NumericInput) -> Scalar:
    r"""Return ``rho^[k](p^a) / p^(a eps)``."""
    require_prime(p)
    require_integer(a, "a", minimum=0)
    require_kappa(k)
    require_finite(eps, "eps")
    F = current_field()
    return rho_pp(p, a, k) / F.number(p) ** (F.number(a) * F.number(eps))


def square_entry_gap_k5(p: int) -> Scalar:
    r"""Return ``rho^[5](p^2) - rho^[5](p)^2``.

    For every prime this is positive by the identity

    ``5((1+x+x^2)^5-(1+x)^5)-((1+x)^5-1)^2``
    ``= x^6(4x^4+15x^3+30x^2+30x+15)``.
    """
    require_prime(p)
    return rho_pp(p, 2, 5) - rho_pp(p, 1, 5) ** 2
