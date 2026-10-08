"""Runtime checks for the domains of the public numerical functions.

Dynamic mpmath values do not relax the mathematical domains: integer arguments
remain integers, prime-power exponents remain nonnegative where required, and
the deformation parameter satisfies ``kappa > 1``.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import mpmath
from sympy import isprime

from ._types import NumericInput


def require_integer(value: Any, name: str, *, minimum: int | None = None) -> int:
    """Return ``value`` as an integer after enforcing an optional lower bound."""
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return int(value)


def require_prime(value: Any, name: str = "p") -> int:
    """Return ``value`` after checking that it is a positive prime integer."""
    prime = require_integer(value, name, minimum=2)
    if not bool(isprime(prime)):
        raise ValueError(f"{name} must be prime")
    return prime


def _numeric_bounds(value: NumericInput, name: str) -> tuple[mpmath.mpf, mpmath.mpf]:
    """Return the represented lower and upper bounds as ``mpf`` values.

    Point values are degenerate intervals; ``mpi`` values contribute their two
    endpoints.
    """
    try:
        if hasattr(value, "a") and hasattr(value, "b"):
            lo = mpmath.mpf(value.a)
            hi = mpmath.mpf(value.b)
        else:
            lo = hi = mpmath.mpf(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a real numeric value") from exc
    if not mpmath.isfinite(lo) or not mpmath.isfinite(hi):
        raise ValueError(f"{name} must be finite")
    return lo, hi


def require_finite(value: NumericInput, name: str) -> None:
    """Require a finite point value or finite interval."""
    _numeric_bounds(value, name)


def require_real_greater_than(value: NumericInput, bound: NumericInput, name: str) -> None:
    """Require every represented value to lie strictly above ``bound``."""
    lo, _ = _numeric_bounds(value, name)
    _, bound_hi = _numeric_bounds(bound, "bound")
    if lo <= bound_hi:
        raise ValueError(f"{name} must be greater than {bound_hi}")


def require_real_at_least(value: NumericInput, bound: NumericInput, name: str) -> None:
    """Require every represented value to lie at or above ``bound``."""
    lo, _ = _numeric_bounds(value, name)
    _, bound_hi = _numeric_bounds(bound, "bound")
    if lo < bound_hi:
        raise ValueError(f"{name} must be at least {bound_hi}")


def require_positive(value: NumericInput, name: str) -> None:
    require_real_greater_than(value, 0, name)


def require_kappa(value: NumericInput, name: str = "kappa") -> None:
    require_real_greater_than(value, 1, name)


def require_epsilon_in_ca_domain(kappa: NumericInput, epsilon: NumericInput) -> None:
    """Require ``epsilon > -1/kappa`` for the reciprocal-CA construction."""
    require_kappa(kappa)
    _, k_hi = _numeric_bounds(kappa, "kappa")
    eps_lo, _ = _numeric_bounds(epsilon, "epsilon")
    boundary = -1 / k_hi
    if eps_lo <= boundary:
        raise ValueError(f"epsilon must be greater than -1/kappa ({boundary})")


def require_precision(value: Any, name: str = "dps", *, minimum: int = 2) -> int:
    return require_integer(value, name, minimum=minimum)


def validated_factorization(mapping: Mapping[int, int]) -> tuple[tuple[int, int], ...]:
    """Validate a prime-factor mapping and return it in increasing prime order.

    Exponents must be positive. The empty mapping represents 1.
    """
    pairs: list[tuple[int, int]] = []
    for raw_prime, raw_exponent in mapping.items():
        prime = require_prime(raw_prime, "factorization prime")
        exponent = require_integer(raw_exponent, "factorization exponent", minimum=1)
        pairs.append((prime, exponent))
    return tuple(sorted(pairs))
