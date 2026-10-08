"""Independently verify a fixed-parameter normalized-endpoint certificate.

The verifier reconstructs the finite Mertens tail, local upper hulls, endpoint
transitions, Robin margins, and initial integer scan from manuscript formulas
and certificate witnesses. It does not import the generator or the package's
theorem-specific arithmetic. Rounded decimal endpoint strings are diagnostic
only; verification recomputes the directed intervals.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from math import isqrt
from pathlib import Path
from typing import Any, cast

import mpmath
from sympy import nextprime, primerange

from .._serialization import read_json
from .runtime import unbounded_decimal_integers

CERTIFICATE_ANALYTIC_BASE = 58_271
TAIL_COEFFICIENT_NUMERATOR = 3375
TAIL_COEFFICIENT_DENOMINATOR = 1022


def _certificate_int(value: Any, name: str, *, minimum: int | None = None) -> int:
    """Read an integer-valued certificate field without permissive coercion."""
    if isinstance(value, bool) or not isinstance(value, int):
        raise AssertionError(f"{name} must be an integer")
    if minimum is not None and value < minimum:
        raise AssertionError(f"{name} must be at least {minimum}")
    return cast(int, value)


def _certificate_decimal_integer(value: Any, name: str, *, minimum: int = 1) -> int:
    """Read a canonical positive decimal-integer string from a certificate."""
    if not isinstance(value, str) or not value or not value.isascii() or not value.isdecimal():
        raise AssertionError(f"{name} must be a decimal integer string")
    if len(value) > 1 and value.startswith("0"):
        raise AssertionError(f"{name} must use canonical decimal notation")
    result = int(value)
    if result < minimum:
        raise AssertionError(f"{name} must be at least {minimum}")
    return result


@contextmanager
def _working_precision(dps: int) -> Iterator[None]:
    """Set both mpmath point and interval precision for this verifier only."""
    if isinstance(dps, bool) or not isinstance(dps, int) or dps < 2:
        raise AssertionError("precision_dps must be an integer at least 2")
    prior_mp, prior_iv = mpmath.mp.dps, mpmath.iv.dps
    mpmath.mp.dps = dps
    mpmath.iv.dps = dps
    try:
        yield
    finally:
        mpmath.mp.dps = prior_mp
        mpmath.iv.dps = prior_iv


def _sigma_minus_one_pp(p: int, a: int) -> Any:
    """Direct interval formula for sigma_{-1}(p^a)."""
    if a < 0:
        raise ValueError("prime-power exponent must be nonnegative")
    p_iv = mpmath.mpi(p)
    one = mpmath.mpi(1)
    return (one - p_iv ** (-(a + 1))) / (one - one / p_iv)


def _rksigma_pp(p: int, a: int, kappa: str) -> Any:
    """Direct interval formula for sigma_{-1}^{[kappa]}(p^a)."""
    if a == 0:
        return mpmath.mpi(1)
    kappa_iv = mpmath.mpi(kappa)
    return _sigma_minus_one_pp(p, a) ** kappa_iv - _sigma_minus_one_pp(p, a - 1) ** kappa_iv


def _rho_pp(p: int, a: int, kappa: str) -> Any:
    """Direct interval formula for the normalized local factor Q_kappa(p^a)."""
    if a == 0:
        return mpmath.mpi(1)
    p_iv = mpmath.mpi(p)
    kappa_iv = mpmath.mpi(kappa)
    one = mpmath.mpi(1)
    base = p_iv**a * _rksigma_pp(p, a, kappa) / kappa_iv
    return base ** (one / (kappa_iv - one))


def _finite_tail_check(kappa: str, claimed: dict[str, Any], dps: int) -> dict[str, Any]:
    iv = mpmath.iv
    kappa_point = mpmath.mpf(kappa)
    a_floor = int(mpmath.floor(3 * (kappa_point + 3) / 4))
    if a_floor != _certificate_int(claimed["A_K_floor"], "finite_tail.A_K_floor", minimum=1):
        raise AssertionError("A_K floor mismatch")

    with _working_precision(dps):
        c_k = mpmath.mpi(0)
        for r in primerange(2, a_floor + 1):
            c_k += iv.log(r)
        x_second = (mpmath.mpi(TAIL_COEFFICIENT_NUMERATOR) / TAIL_COEFFICIENT_DENOMINATOR) * (
            c_k + mpmath.mpi(1) / 3
        )
        analytic_cutoff = max(CERTIFICATE_ANALYTIC_BASE, int(mpmath.ceil(mpmath.mpf(x_second.b))))
        if analytic_cutoff != _certificate_int(
            claimed["analytic_cutoff"], "finite_tail.analytic_cutoff", minimum=2
        ):
            raise AssertionError("analytic cutoff mismatch")

        primes = list(primerange(2, analytic_cutoff + 1))
        theta = mpmath.mpi(0)
        theta_third = mpmath.mpi(0)
        log_mertens = mpmath.mpi(0)
        third_index = 0
        last_bad: int | None = None
        terminal_min: mpmath.mpf | None = None
        terminal_min_prime: int | None = None

        for p in primes:
            p_iv = mpmath.mpi(p)
            theta += iv.log(p_iv)
            log_mertens -= iv.log(1 - 1 / p_iv)
            bound = (p - 1) / 3
            while third_index < len(primes) and primes[third_index] <= bound:
                theta_third += iv.log(primes[third_index])
                third_index += 1
            log_n_lower = theta + theta_third - c_k
            if log_n_lower.a <= 1:
                last_bad = p
                terminal_min = None
                terminal_min_prime = None
                continue
            margin = iv.euler + iv.log(iv.log(log_n_lower)) - log_mertens
            if margin.a <= 0:
                last_bad = p
                terminal_min = None
                terminal_min_prime = None
                continue
            lower = mpmath.mpf(margin.a)
            if terminal_min is None or lower < terminal_min:
                terminal_min = lower
                terminal_min_prime = p

        expected_tail = 2 if last_bad is None else int(nextprime(last_bad))
        if expected_tail > analytic_cutoff:
            expected_tail = int(nextprime(analytic_cutoff))
        if last_bad != claimed["last_bad_prime"]:
            raise AssertionError("last bad finite-tail prime mismatch")
        if expected_tail != _certificate_int(
            claimed["tail_start_prime"], "finite_tail.tail_start_prime", minimum=2
        ):
            raise AssertionError("finite-tail start prime mismatch")
        if len(primes) != _certificate_int(
            claimed["prime_rows_checked"], "finite_tail.prime_rows_checked", minimum=1
        ):
            raise AssertionError("finite-tail prime count mismatch")
        if expected_tail <= analytic_cutoff:
            if terminal_min is None or terminal_min <= 0:
                raise AssertionError("terminal finite-tail block is not certified positive")
            if terminal_min_prime != claimed["terminal_block_min_margin_prime"]:
                raise AssertionError("terminal finite-tail minimum prime mismatch")

        return {
            "tail_start_prime": expected_tail,
            "last_bad_prime": last_bad,
            "terminal_min_margin": terminal_min,
            "terminal_min_prime": terminal_min_prime,
            "prime_rows_checked": len(primes),
        }


def _entry_check(kappa: str, claimed: dict[str, Any], dps: int) -> Any:
    iv = mpmath.iv
    q = _certificate_int(claimed["prime"], "entry_parameter.prime", minimum=2)
    exponent = _certificate_int(claimed["exponent"], "entry_parameter.exponent", minimum=1)
    limit = _certificate_int(
        claimed["candidate_limit"], "entry_parameter.candidate_limit", minimum=1
    )
    if exponent > limit:
        raise AssertionError("entry_parameter.exponent exceeds candidate_limit")
    with _working_precision(dps):
        log_q = iv.log(q)
        values = [iv.log(_rho_pp(q, a, kappa)) / (a * log_q) for a in range(1, limit + 1)]
        best = values[exponent - 1]
        for a, value in enumerate(values, start=1):
            if a != exponent and value.b >= best.a:
                raise AssertionError("entry exponent is not uniquely certified")
        tail_upper = iv.log(mpmath.mpi(q) / (q - 1)) / ((limit + 1) * log_q)
        if tail_upper.b >= best.a:
            raise AssertionError("entry exponent tail is not excluded")
        return best


def _verify_local_profiles(
    kappa: str,
    profiles: list[dict[str, Any]],
    epsilon: Any,
    q: int,
    boundary_exponent: int,
    dps: int,
) -> list[dict[str, Any]]:
    iv = mpmath.iv
    expected_primes = list(primerange(2, q + 1))
    claimed_primes = [
        _certificate_int(row["prime"], "local_profiles[].prime", minimum=2) for row in profiles
    ]
    if claimed_primes != expected_primes:
        raise AssertionError("local-profile prime list is incomplete")

    reconstructed: list[dict[str, Any]] = []
    with _working_precision(dps):
        for row in profiles:
            p = _certificate_int(row["prime"], "local_profiles[].prime", minimum=2)
            cap = _certificate_int(
                row["exponent_cap"], f"local profile p={p} exponent_cap", minimum=1
            )
            raw_hull = row["hull_vertices"]
            if not isinstance(raw_hull, list):
                raise AssertionError(f"hull vertices must be a list at p={p}")
            hull = [_certificate_int(a, f"hull vertex at p={p}", minimum=0) for a in raw_hull]
            if not hull or hull[0] != 0 or hull[-1] != cap:
                raise AssertionError(f"invalid hull endpoints at p={p}")
            if any(left >= right for left, right in zip(hull, hull[1:])):
                raise AssertionError(f"hull vertices are not strictly increasing at p={p}")

            # For a >= cap, Q(p^(a+1))/Q(p^a) is bounded by
            # S_(a+1)/S_(a-1), which decreases with a.
            ratio_bound = _sigma_minus_one_pp(p, cap + 1) / _sigma_minus_one_pp(p, cap - 1)
            tail_slope = iv.log(ratio_bound) / iv.log(p)
            if tail_slope.b >= epsilon.a:
                raise AssertionError(f"local exponent tail not excluded at p={p}")

            log_p = iv.log(p)
            y = [mpmath.mpi(0)] + [iv.log(_rho_pp(p, a, kappa)) / log_p for a in range(1, cap + 1)]
            slopes: list[Any] = []
            for left, right in zip(hull, hull[1:]):
                slopes.append((y[right] - y[left]) / (right - left))
            for first, second in zip(slopes, slopes[1:]):
                if first.a <= second.b:
                    raise AssertionError(f"hull slopes are not strictly decreasing at p={p}")

            for left, right in zip(hull, hull[1:]):
                for a in range(left + 1, right):
                    chord = y[left] + (a - left) * (y[right] - y[left]) / (right - left)
                    if y[a].b >= chord.a:
                        raise AssertionError(f"point a={a} is not certified below hull at p={p}")

            local_count = 0
            for (left, right), event in zip(zip(hull, hull[1:]), slopes):
                if event.a > epsilon.b:
                    reconstructed.append(
                        {
                            "prime": p,
                            "from_exponent": left,
                            "to_exponent": right,
                            "interval": event,
                        }
                    )
                    local_count += 1
                elif event.b < epsilon.a:
                    continue
                elif p == q and left == 0 and right == boundary_exponent:
                    reconstructed.append(
                        {
                            "prime": p,
                            "from_exponent": left,
                            "to_exponent": right,
                            "interval": event,
                            "boundary_event": True,
                        }
                    )
                    local_count += 1
                else:
                    raise AssertionError(f"event is not separated from boundary epsilon at p={p}")
            if local_count != _certificate_int(
                row["relevant_event_count"], f"local profile p={p} relevant_event_count", minimum=0
            ):
                raise AssertionError(f"relevant event count mismatch at p={p}")

    reconstructed.sort(key=lambda event: mpmath.mpf(event["interval"].a), reverse=True)
    for current, following in zip(reconstructed, reconstructed[1:]):
        if current["interval"].a <= following["interval"].b:
            raise AssertionError("global endpoint event intervals overlap")
    return reconstructed


def _verify_events_and_endpoints(
    kappa: str,
    reconstructed: list[dict[str, Any]],
    claimed_events: list[dict[str, Any]],
    endpoint_data: dict[str, Any],
    dps: int,
) -> dict[str, Any]:
    iv = mpmath.iv
    if len(reconstructed) != len(claimed_events):
        raise AssertionError("global event count mismatch")
    for actual, claimed in zip(reconstructed, claimed_events):
        event_key = (actual["prime"], actual["from_exponent"], actual["to_exponent"])
        claimed_key = (
            _certificate_int(claimed["prime"], "global_events[].prime", minimum=2),
            _certificate_int(claimed["from_exponent"], "global_events[].from_exponent", minimum=0),
            _certificate_int(claimed["to_exponent"], "global_events[].to_exponent", minimum=1),
        )
        if event_key != claimed_key:
            raise AssertionError("global event order mismatch")
        if claimed.get("boundary_event", False) is not actual.get("boundary_event", False):
            raise AssertionError("global boundary-event marker mismatch")

    rows = endpoint_data["rows"]
    if len(rows) != len(reconstructed):
        raise AssertionError("endpoint-row count mismatch")

    factorization: dict[int, int] = {}
    n = 1
    first_above: int | None = None
    min_margin: mpmath.mpf | None = None
    min_n: int | None = None
    log_prime_cache: dict[int, Any] = {}
    log_rho_cache: dict[tuple[int, int], Any] = {}
    previous_n = 0

    with _working_precision(dps):
        for index, (event, row) in enumerate(zip(reconstructed, rows), start=1):
            p = int(event["prime"])
            old = int(event["from_exponent"])
            new = int(event["to_exponent"])
            if factorization.get(p, 0) != old:
                raise AssertionError(f"endpoint state mismatch at event {index}")
            n *= p ** (new - old)
            factorization[p] = new
            if n <= previous_n:
                raise AssertionError(f"endpoint sequence is not increasing at event {index}")
            previous_n = n

            row_key = (
                _certificate_int(
                    row["event_prime"], f"endpoint row {index} event_prime", minimum=2
                ),
                _certificate_int(
                    row["from_exponent"], f"endpoint row {index} from_exponent", minimum=0
                ),
                _certificate_int(
                    row["to_exponent"], f"endpoint row {index} to_exponent", minimum=1
                ),
            )
            if row_key != (p, old, new):
                raise AssertionError(f"endpoint event metadata mismatch at event {index}")
            row_index = _certificate_int(row["index"], f"endpoint row {index} index", minimum=1)
            row_endpoint = _certificate_decimal_integer(
                row["endpoint"], f"endpoint row {index} endpoint"
            )
            if row_index != index or row_endpoint != n:
                raise AssertionError(f"endpoint witness mismatch at event {index}")
            if _certificate_int(
                row["largest_prime"], f"endpoint row {index} largest_prime", minimum=2
            ) != max(factorization):
                raise AssertionError(f"largest-prime witness mismatch at event {index}")

            if n > 5040:
                log_n = mpmath.mpi(0)
                log_q = mpmath.mpi(0)
                for prime, exponent in factorization.items():
                    log_prime = log_prime_cache.get(prime)
                    if log_prime is None:
                        log_prime = iv.log(prime)
                        log_prime_cache[prime] = log_prime
                    log_n += exponent * log_prime

                    local_factor_key = (prime, exponent)
                    log_rho = log_rho_cache.get(local_factor_key)
                    if log_rho is None:
                        log_rho = iv.log(_rho_pp(prime, exponent, kappa))
                        log_rho_cache[local_factor_key] = log_rho
                    log_q += log_rho
                margin = iv.euler + iv.log(iv.log(log_n)) - log_q
                if margin.a <= 0:
                    raise AssertionError(f"Robin endpoint margin failed at N={n}")
                lower = mpmath.mpf(margin.a)
                if first_above is None:
                    first_above = n
                if min_margin is None or lower < min_margin:
                    min_margin = lower
                    min_n = n

    if first_above is None:
        raise AssertionError("no endpoint above 5040")
    if endpoint_data.get("all_positive") is not True:
        raise AssertionError("endpoint all-positive summary is not true")
    if _certificate_int(
        endpoint_data["event_count"], "endpoint_verification.event_count", minimum=1
    ) != len(reconstructed):
        raise AssertionError("endpoint event-count summary mismatch")
    if (
        _certificate_decimal_integer(
            endpoint_data["first_endpoint_above_5040"],
            "endpoint_verification.first_endpoint_above_5040",
        )
        != first_above
    ):
        raise AssertionError("first endpoint above 5040 mismatch")
    if (
        _certificate_decimal_integer(
            endpoint_data["tail_endpoint"], "endpoint_verification.tail_endpoint"
        )
        != n
    ):
        raise AssertionError("tail endpoint mismatch")
    if (
        _certificate_decimal_integer(
            endpoint_data["minimum_log_margin_endpoint"],
            "endpoint_verification.minimum_log_margin_endpoint",
        )
        != min_n
    ):
        raise AssertionError("minimum endpoint-margin location mismatch")
    return {
        "event_count": len(reconstructed),
        "first_endpoint_above_5040": first_above,
        "tail_endpoint": n,
        "minimum_log_margin": min_margin,
        "minimum_log_margin_endpoint": min_n,
    }


def _smallest_prime_factor_table(limit: int) -> list[int]:
    """Return smallest-prime-factor data for every integer through ``limit``."""
    if limit < 1:
        return [0] * (limit + 1)
    spf = list(range(limit + 1))
    spf[0] = 0
    spf[1] = 1
    for p in range(2, isqrt(limit) + 1):
        if spf[p] != p:
            continue
        for multiple in range(p * p, limit + 1, p):
            if spf[multiple] == multiple:
                spf[multiple] = p
    return spf


def _verify_initial_scan(
    kappa: str,
    claimed: dict[str, Any],
    first_endpoint: int,
    dps: int,
) -> dict[str, Any]:
    iv = mpmath.iv
    lower = 5041
    upper = first_endpoint - 1
    claimed_lower = _certificate_int(claimed["lower"], "initial_integer_scan.lower", minimum=1)
    claimed_upper = _certificate_int(claimed["upper"], "initial_integer_scan.upper", minimum=0)
    if claimed_lower != lower or claimed_upper != upper:
        raise AssertionError("initial scan range mismatch")
    if _certificate_int(claimed["count"], "initial_integer_scan.count", minimum=0) != max(
        0, upper - lower + 1
    ):
        raise AssertionError("initial scan count mismatch")
    if claimed.get("all_positive") is not True:
        raise AssertionError("initial-scan all-positive summary is not true")

    minimum: mpmath.mpf | None = None
    minimum_n: int | None = None
    spf = _smallest_prime_factor_table(upper)
    with _working_precision(dps):
        e_gamma = iv.exp(iv.euler)
        kappa_iv = mpmath.mpi(kappa)
        one = mpmath.mpi(1)
        root = one / (kappa_iv - one)
        local_base_cache: dict[tuple[int, int], Any] = {}
        for n in range(lower, upper + 1):
            remainder = n
            base_product = one
            while remainder > 1:
                p = spf[remainder]
                a = 0
                while remainder % p == 0:
                    remainder //= p
                    a += 1
                key = (p, a)
                local_base = local_base_cache.get(key)
                if local_base is None:
                    p_iv = mpmath.mpi(p)
                    local_base = p_iv**a * _rksigma_pp(p, a, kappa) / kappa_iv
                    local_base_cache[key] = local_base
                base_product *= local_base
            margin = e_gamma * iv.log(iv.log(n)) - base_product**root
            if margin.a <= 0:
                raise AssertionError(f"initial integer scan failed at n={n}")
            value = mpmath.mpf(margin.a)
            if minimum is None or value < minimum:
                minimum = value
                minimum_n = n
    if minimum_n is None:
        if claimed["minimum_at_n"] is not None:
            raise AssertionError("empty initial scan must not claim a minimum location")
    elif minimum_n != _certificate_int(
        claimed["minimum_at_n"], "initial_integer_scan.minimum_at_n", minimum=lower
    ):
        raise AssertionError("initial scan minimum location mismatch")
    return {
        "count": max(0, upper - lower + 1),
        "minimum_margin": minimum,
        "minimum_at_n": minimum_n,
    }


def verify_certificate(certificate: Path) -> dict[str, Any]:
    """Verify one certificate by independent reconstruction.

    A failed check raises ``AssertionError``; success returns a compact report.
    """
    with unbounded_decimal_integers():
        data = read_json(certificate)
        if not isinstance(data, dict):
            raise AssertionError("bounded-parameter certificate must be a JSON object")
        if data.get("schema_version") != "bounded-parameter-endpoints-1.0.0":
            raise AssertionError("unsupported bounded-parameter certificate schema")
        if data.get("status") != "directed-interval-finite-certificate":
            raise AssertionError("unexpected bounded-parameter certificate status")
        dps = _certificate_int(data["precision_dps"], "precision_dps", minimum=2)
        if not isinstance(data.get("kappa"), str):
            raise AssertionError("kappa must be stored as a decimal string")
        kappa = data["kappa"]

        with _working_precision(dps):
            try:
                kappa_point = mpmath.mpf(kappa)
            except (TypeError, ValueError) as exc:
                raise AssertionError("kappa is not a valid real number") from exc
            if not mpmath.isfinite(kappa_point) or kappa_point <= 1:
                raise AssertionError("kappa must be a finite real number greater than 1")
            expected_scope = f"all n > 5040 and all real 1 < kappa <= {kappa}"
            if data.get("theorem_scope") != expected_scope:
                raise AssertionError("theorem scope does not match kappa")

            tail = _finite_tail_check(kappa, data["finite_tail"], dps)
            q = int(tail["tail_start_prime"])
            if q != _certificate_int(
                data["entry_parameter"]["prime"], "entry_parameter.prime", minimum=2
            ):
                raise AssertionError("entry prime does not equal tail-start prime")
            epsilon = _entry_check(kappa, data["entry_parameter"], dps)
            boundary_exponent = _certificate_int(
                data["entry_parameter"]["exponent"],
                "entry_parameter.exponent",
                minimum=1,
            )
            events = _verify_local_profiles(
                kappa,
                data["local_profiles"],
                epsilon,
                q,
                boundary_exponent,
                dps,
            )
            if not events:
                raise AssertionError("certificate contains no canonical endpoint events")
            final_event = events[-1]
            if (
                final_event["prime"],
                final_event["from_exponent"],
                final_event["to_exponent"],
                final_event.get("boundary_event"),
            ) != (q, 0, boundary_exponent, True):
                raise AssertionError("final event is not the certified tail-prime entry event")
            endpoints = _verify_events_and_endpoints(
                kappa,
                events,
                data["global_events"],
                data["endpoint_verification"],
                dps,
            )
            initial = _verify_initial_scan(
                kappa,
                data["initial_integer_scan"],
                int(endpoints["first_endpoint_above_5040"]),
                dps,
            )

    report = {
        "certificate": str(certificate),
        "kappa": kappa,
        "finite_tail_start_prime": q,
        "finite_tail_last_bad_prime": tail["last_bad_prime"],
        "finite_tail_prime_rows": tail["prime_rows_checked"],
        "canonical_events": endpoints["event_count"],
        "first_endpoint_above_5040": endpoints["first_endpoint_above_5040"],
        "minimum_endpoint_log_margin_at": endpoints["minimum_log_margin_endpoint"],
        "initial_integer_checks": initial["count"],
        "initial_minimum_at": initial["minimum_at_n"],
        "all_checks_pass": True,
    }
    return report


def main(
    argv: Sequence[str] | None = None,
    *,
    default_certificate: Path | None = None,
) -> int:
    """Verify one bounded-parameter certificate from the command line."""
    parser = argparse.ArgumentParser(
        description="Independently verify a bounded-parameter certificate."
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
    report = verify_certificate(args.certificate)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
