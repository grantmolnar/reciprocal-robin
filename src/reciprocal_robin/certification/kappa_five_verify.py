"""Independently verify the finite certificate for ``1 < kappa <= 5``.

This module reconstructs the manuscript's prime-power formulas directly with
directed interval arithmetic. It does not import the generator or the package's
theorem-specific arithmetic. Decimal endpoint strings are diagnostic only; the
verifier makes its decisions from freshly computed intervals and then compares
the reconstructed object with the checked-in JSON.
"""

from __future__ import annotations

import argparse
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from fractions import Fraction
from pathlib import Path
from typing import Any, cast

import mpmath
from sympy import primerange

from .._serialization import json_text, read_json


@contextmanager
def _working_precision(dps: int) -> Iterator[None]:
    old_mp = mpmath.mp.dps
    old_iv = mpmath.iv.dps
    mpmath.mp.dps = dps
    mpmath.iv.dps = dps
    try:
        yield
    finally:
        mpmath.mp.dps = old_mp
        mpmath.iv.dps = old_iv


def _sigma_minus_one_pp(p: int, a: int) -> Any:
    p_iv = mpmath.mpi(p)
    one = mpmath.mpi(1)
    return (one - p_iv ** (-(a + 1))) / (one - one / p_iv)


def _rksigma_pp(p: int, a: int, kappa: int) -> Any:
    if a == 0:
        return mpmath.mpi(1)
    return _sigma_minus_one_pp(p, a) ** kappa - _sigma_minus_one_pp(p, a - 1) ** kappa


def _rho_pp(p: int, a: int, kappa: int) -> Any:
    if a == 0:
        return mpmath.mpi(1)
    p_iv = mpmath.mpi(p)
    kappa_iv = mpmath.mpi(kappa)
    base = p_iv**a * _rksigma_pp(p, a, kappa) / kappa_iv
    return base ** (mpmath.mpi(1) / (kappa_iv - 1))


def _normalized_local_objective(p: int, a: int, kappa: int, eps: Any) -> Any:
    return _rho_pp(p, a, kappa) / mpmath.mpi(p) ** (a * eps)


MERTENS_EXTENSION_LOWER = 19_421
DUSART_ANALYTIC_THRESHOLD = 2_278_382
MERTENS_CORRECTION = "0.10952"
SCHEMA_VERSION = "kappa-five-interval-1.1.0"
CERTIFICATE_STATUS = "directed-interval-finite-certificate"
THEOREM_SCOPE = (
    "shared finite inputs for the unconditional interval 1 < kappa <= 5 "
    "and the uniform explicit upper bound"
)


def _certificate_int(value: Any, name: str, *, minimum: int | None = None) -> int:
    """Read an integer-valued certificate field without permissive coercion."""
    if isinstance(value, bool) or not isinstance(value, int):
        raise AssertionError(f"{name} must be an integer")
    if minimum is not None and value < minimum:
        raise AssertionError(f"{name} must be at least {minimum}")
    return cast(int, value)


def _lo(x: Any, digits: int) -> str:
    return str(mpmath.nstr(mpmath.mpf(x.a), digits))


def _hi(x: Any, digits: int) -> str:
    return str(mpmath.nstr(mpmath.mpf(x.b), digits))


def _interval_fraction(value: Fraction) -> Any:
    return mpmath.mpi(value.numerator) / mpmath.mpi(value.denominator)


def _reconstruct_certificate(dps: int = 60) -> dict[str, Any]:
    """Independently reconstruct the full finite certificate at ``dps`` digits."""
    iv = mpmath.iv
    with _working_precision(dps):
        e_gamma = iv.exp(iv.euler)

        theta = mpmath.mpi(0)
        mertens = Fraction(1, 1)
        rows: list[dict[str, Any]] = []
        min_margin = None
        min_prime = None
        for p in primerange(2, MERTENS_EXTENSION_LOWER + 1):
            theta += iv.log(p)
            if p >= MERTENS_EXTENSION_LOWER:
                continue
            mertens *= Fraction(p, p - 1)
            if p < 11:
                continue
            margin = e_gamma * iv.log(2 * theta) - _interval_fraction(mertens)
            certified = margin.a > 0
            rows.append(
                {
                    "P": int(p),
                    "margin_lo": _lo(margin, dps),
                    "margin_hi": _hi(margin, dps),
                    "certified_positive": bool(certified),
                }
            )
            lower = mpmath.mpf(margin.a)
            if min_margin is None or lower < min_margin:
                min_margin = lower
                min_prime = int(p)

        theta_at_threshold = theta
        primorial_threshold_margin = mpmath.mpi(19_183) - theta_at_threshold

        # The Mertens product is constant between primes, while the right-hand
        # side increases here. It is enough to check prime endpoints. Summing
        # logarithms avoids an enormous exact rational product.
        log_mertens = mpmath.mpi(0)
        correction = mpmath.mpi(MERTENS_CORRECTION)
        extension_count = 0
        extension_all_verified = True
        extension_min_margin = None
        extension_min_prime = None
        extension_max_required = None
        extension_max_prime = None
        for p in primerange(2, DUSART_ANALYTIC_THRESHOLD):
            p_iv = mpmath.mpi(p)
            log_mertens -= iv.log(1 - 1 / p_iv)
            if p < MERTENS_EXTENSION_LOWER:
                continue
            extension_count += 1
            L = iv.log(p_iv)
            rhs_log = iv.euler + iv.log(L) + iv.log(1 + correction / L**2)
            log_margin = rhs_log - log_mertens
            certified = bool(log_margin.a > 0)
            extension_all_verified = extension_all_verified and certified
            lower_margin = mpmath.mpf(log_margin.a)
            if extension_min_margin is None or lower_margin < extension_min_margin:
                extension_min_margin = lower_margin
                extension_min_prime = int(p)

            required = L**2 * (iv.exp(log_mertens - iv.euler) / L - 1)
            required_upper = mpmath.mpf(required.b)
            if extension_max_required is None or required_upper > extension_max_required:
                extension_max_required = required_upper
                extension_max_prime = int(p)

        eps0 = iv.log(_rho_pp(3, 3, 5) / _rho_pp(3, 2, 5)) / iv.log(3)

        def obj(p: int, a: int) -> Any:
            return _normalized_local_objective(p, a, 5, eps0)

        obj2 = {a: obj(2, a) for a in range(0, 6)}
        obj3 = {a: obj(3, a) for a in range(0, 5)}
        obj5 = {a: obj(5, a) for a in range(0, 4)}

        tail2 = mpmath.mpi(2) / (mpmath.mpi(2) ** (6 * eps0))
        tail3 = mpmath.mpi(3) / mpmath.mpi(2) / (mpmath.mpi(3) ** (5 * eps0))
        tail5 = mpmath.mpi(5) / mpmath.mpi(4) / (mpmath.mpi(5) ** (3 * eps0))
        excluded7 = mpmath.mpi(7) / mpmath.mpi(6) / (mpmath.mpi(7) ** (2 * eps0))

        seed_checks = {
            "p2_a4_beats_a0_to_a5_except_a4": all(obj2[4].a > obj2[a].b for a in obj2 if a != 4),
            "p2_a4_beats_tail_a_ge_6": obj2[4].a > tail2.b,
            # The a=2 and a=3 values are exactly equal by the definition of
            # eps0. Interval overlap checks the numerical reconstruction only.
            "p3_a2_a3_tie_exact_by_definition": True,
            "p3_a2_a3_interval_consistency": not (obj3[2].b < obj3[3].a or obj3[3].b < obj3[2].a),
            "p3_a3_beats_a0_a1_a4": all(obj3[3].a > obj3[a].b for a in (0, 1, 4)),
            "p3_a3_beats_tail_a_ge_5": obj3[3].a > tail3.b,
            "p5_a2_beats_a0_a1_a3": all(obj5[2].a > obj5[a].b for a in (0, 1, 3)),
            "p5_a2_beats_tail_a_ge_3": obj5[2].a > tail5.b,
            "all_p_ge_7_excluded": excluded7.b < 1,
        }

        low_support_margin = e_gamma * iv.log(20 * iv.log(2)) - mpmath.mpi(35) / 8

        seed_factorization = ((2, 4), (3, 3), (5, 2))
        seed_endpoint = 1
        for p, exponent in seed_factorization:
            seed_endpoint *= p**exponent

        return {
            "schema_version": SCHEMA_VERSION,
            "status": CERTIFICATE_STATUS,
            "precision_dps": dps,
            "theorem_scope": THEOREM_SCOPE,
            "mertens_margin": {
                "range": "prime P with 11 <= P < 19421",
                "all_positive": all(row["certified_positive"] for row in rows),
                "minimum_lower_margin": mpmath.nstr(min_margin, dps),
                "minimum_at_prime": min_prime,
                "rows": rows,
            },
            "primorial_threshold": {
                "quantity": "19183 - theta(19421)",
                "theta_upper": _hi(theta_at_threshold, dps),
                "margin_lo": _lo(primorial_threshold_margin, dps),
                "margin_hi": _hi(primorial_threshold_margin, dps),
                "certified_positive": bool(primorial_threshold_margin.a > 0),
            },
            "mertens_extension": {
                "range": (
                    "prime P with 19421 <= P < 2278382; real x follows "
                    "by monotonicity between prime endpoints"
                ),
                "correction_constant": MERTENS_CORRECTION,
                "prime_count": extension_count,
                "all_verified": bool(extension_all_verified),
                "minimum_log_margin_lower": mpmath.nstr(extension_min_margin, dps),
                "minimum_log_margin_at_prime": extension_min_prime,
                "maximum_required_correction_upper": mpmath.nstr(extension_max_required, dps),
                "maximum_required_correction_at_prime": extension_max_prime,
            },
            "normalized_ca_seed": {
                "epsilon_lo": _lo(eps0, dps),
                "epsilon_hi": _hi(eps0, dps),
                "largest_endpoint": seed_endpoint,
                "factorization": {str(p): exponent for p, exponent in seed_factorization},
                "exact_tie_identity": (
                    "Q(3^3)*3^(-3*epsilon0) = Q(3^2)*3^(-2*epsilon0) "
                    "by the definition of epsilon0"
                ),
                "checks": seed_checks,
                "all_checks_pass": all(seed_checks.values()),
            },
            "low_support_margin": {
                "quantity": "exp(gamma)*log(20*log(2)) - 35/8",
                "margin_lo": _lo(low_support_margin, dps),
                "margin_hi": _hi(low_support_margin, dps),
                "certified_positive": bool(low_support_margin.a > 0),
            },
            # The manuscript proves this identity algebraically; the string is
            # recorded as provenance, not as an interval check.
            "symbolic_square_entry_identity": ("x^6*(4*x^4+15*x^3+30*x^2+30*x+15) > 0"),
        }


def verify_certificate(certificate: Path) -> dict[str, Any]:
    """Independently reconstruct and verify one ``kappa <= 5`` certificate."""
    checked_in = read_json(certificate)
    if not isinstance(checked_in, dict):
        raise AssertionError("kappa<=5 certificate must be a JSON object")
    if checked_in.get("schema_version") != SCHEMA_VERSION:
        raise AssertionError("unsupported kappa<=5 certificate schema")
    if checked_in.get("status") != CERTIFICATE_STATUS:
        raise AssertionError("unexpected kappa<=5 certificate status")
    if checked_in.get("theorem_scope") != THEOREM_SCOPE:
        raise AssertionError("unexpected kappa<=5 theorem scope")
    dps = _certificate_int(checked_in.get("precision_dps"), "precision_dps", minimum=2)

    regenerated = _reconstruct_certificate(dps=dps)
    matches = checked_in == regenerated
    finite_checks_pass = (
        regenerated["mertens_margin"]["all_positive"]
        and regenerated["primorial_threshold"]["certified_positive"]
        and regenerated["mertens_extension"]["all_verified"]
        and regenerated["normalized_ca_seed"]["all_checks_pass"]
        and regenerated["low_support_margin"]["certified_positive"]
    )
    if not finite_checks_pass:
        raise AssertionError("independent finite checks failed")
    if not matches:
        raise AssertionError("checked-in kappa<=5 artifact differs from reconstruction")
    return {
        "certificate": str(certificate),
        "finite_checks_pass": True,
        "mertens_rows": len(regenerated["mertens_margin"]["rows"]),
        "minimum_margin_prime": regenerated["mertens_margin"]["minimum_at_prime"],
        "primorial_threshold_pass": regenerated["primorial_threshold"]["certified_positive"],
        "mertens_extension_primes": regenerated["mertens_extension"]["prime_count"],
        "maximum_required_correction_prime": regenerated["mertens_extension"][
            "maximum_required_correction_at_prime"
        ],
        "normalized_endpoint": regenerated["normalized_ca_seed"]["largest_endpoint"],
        "matches_checked_in_artifact": True,
    }


def main(
    argv: Sequence[str] | None = None,
    *,
    default_certificate: Path | None = None,
) -> int:
    parser = argparse.ArgumentParser(
        description="Independently verify the directed-interval kappa<=5 certificate."
    )
    parser.add_argument(
        "--certificate",
        type=Path,
        default=default_certificate,
        required=default_certificate is None,
        help=(
            "certificate JSON path"
            + (f" (default: {default_certificate})" if default_certificate is not None else "")
        ),
    )
    args = parser.parse_args(argv)
    print(json_text(verify_certificate(args.certificate)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
