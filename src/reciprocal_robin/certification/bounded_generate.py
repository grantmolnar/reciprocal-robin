"""Generate a directed-interval certificate at a fixed reciprocal parameter.

The certificate combines the analytic tail theorem from the manuscript with four
finite checks: a repeated-support/Mertens tail handoff, certified local hulls and
endpoint transitions, direct verification of the remaining endpoints above 5040,
and a direct scan from 5041 to the first such endpoint. Monotonicity in ``kappa``
then extends a certificate at ``K`` to ``1 < kappa <= K``.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from dataclasses import dataclass
from math import isqrt
from pathlib import Path
from typing import Any

import mpmath
from sympy import nextprime, primerange

from .. import MpfField, MpiField, working_precision, using_field
from .._serialization import write_json
from ..arithmetic import rksigma_pp, sigma_minus_one_pp
from ..normalized_ca import rho_pp
from .runtime import unbounded_decimal_integers

CERTIFICATE_ANALYTIC_BASE = 58_271
TAIL_COEFFICIENT_NUMERATOR = 3375
TAIL_COEFFICIENT_DENOMINATOR = 1022


def _lo(x: Any, digits: int) -> str:
    """Render the lower endpoint for diagnostics, not as a decimal enclosure."""
    return str(mpmath.nstr(mpmath.mpf(x.a), digits))


def _hi(x: Any, digits: int) -> str:
    """Render the upper endpoint for diagnostics, not as a decimal enclosure."""
    return str(mpmath.nstr(mpmath.mpf(x.b), digits))


def _point(x: Any) -> mpmath.mpf:
    """Return a point representative used only to propose a candidate."""
    if hasattr(x, "a") and hasattr(x, "b"):
        return (mpmath.mpf(x.a) + mpmath.mpf(x.b)) / 2
    return mpmath.mpf(x)


def certificate_filename(kappa: str) -> str:
    """Return the conventional filename for a fixed-parameter certificate."""
    safe = kappa.replace(".", "p").replace("-", "m")
    return f"kappa_{safe}_endpoint_certificate.json"


def default_output_path(kappa: str, directory: Path | None = None) -> Path:
    """Return the default output path used by the installed CLI."""
    base = Path.cwd() if directory is None else directory
    return base / certificate_filename(kappa)


def _theta_cutoff(kappa: mpmath.mpf) -> int:
    """Return the integral repeated-support cutoff used in the tail bound."""
    return int(mpmath.floor(3 * (kappa + 3) / 4))


@dataclass(frozen=True)
class _CertifiedEvent:
    """One interval-certified transition in the canonical endpoint sequence."""

    prime: int
    from_exponent: int
    to_exponent: int
    interval: Any
    boundary_event: bool = False

    def record(self, dps: int) -> dict[str, Any]:
        """Return the stable JSON record for this event."""
        row: dict[str, Any] = {
            "prime": self.prime,
            "from_exponent": self.from_exponent,
            "to_exponent": self.to_exponent,
            "epsilon_lo": _lo(self.interval, dps),
            "epsilon_hi": _hi(self.interval, dps),
        }
        if self.boundary_event:
            row["boundary_event"] = True
        return row


def _finite_tail_scan(kappa_text: str, dps: int) -> dict[str, Any]:
    """Certify the repeated-support/Mertens tail and return its first safe prime."""
    iv = mpmath.iv
    kappa_point = mpmath.mpf(kappa_text)
    a_floor = _theta_cutoff(kappa_point)

    with working_precision(dps), using_field(MpiField()):
        c_k = mpmath.mpi(0)
        for r in primerange(2, a_floor + 1):
            c_k += iv.log(r)

        # For P >= 58271, the manuscript reduces the fixed-K tail to
        # (1022/3375) P > C_K + 1/3.
        x_second = (mpmath.mpi(TAIL_COEFFICIENT_NUMERATOR) / TAIL_COEFFICIENT_DENOMINATOR) * (
            c_k + mpmath.mpi(1) / 3
        )
        analytic_cutoff = max(CERTIFICATE_ANALYTIC_BASE, int(mpmath.ceil(mpmath.mpf(x_second.b))))
        primes = list(primerange(2, analytic_cutoff + 1))

        theta = mpmath.mpi(0)
        theta_third = mpmath.mpi(0)
        log_mertens = mpmath.mpi(0)
        third_index = 0
        last_bad: int | None = None
        post_bad_min: mpmath.mpf | None = None
        post_bad_min_prime: int | None = None
        rows_checked = 0

        for p in primes:
            p_iv = mpmath.mpi(p)
            theta += iv.log(p_iv)
            log_mertens -= iv.log(1 - 1 / p_iv)
            third_bound = (p - 1) / 3
            while third_index < len(primes) and primes[third_index] <= third_bound:
                theta_third += iv.log(primes[third_index])
                third_index += 1

            log_n_lower = theta + theta_third - c_k
            rows_checked += 1
            if log_n_lower.a <= 1:
                last_bad = p
                post_bad_min = None
                post_bad_min_prime = None
                continue

            margin = iv.euler + iv.log(iv.log(log_n_lower)) - log_mertens
            if margin.a <= 0:
                last_bad = p
                post_bad_min = None
                post_bad_min_prime = None
                continue

            lower = mpmath.mpf(margin.a)
            if post_bad_min is None or lower < post_bad_min:
                post_bad_min = lower
                post_bad_min_prime = p

        if last_bad is None:
            tail_start = 2
        else:
            tail_start = int(nextprime(last_bad))
        if tail_start > analytic_cutoff:
            # No positive terminal block occurs before the analytic cutoff.
            tail_start = int(nextprime(analytic_cutoff))
            post_bad_min = None
            post_bad_min_prime = None

        return {
            "A_K_floor": a_floor,
            "C_K_lo": _lo(c_k, dps),
            "C_K_hi": _hi(c_k, dps),
            "analytic_cutoff": analytic_cutoff,
            "prime_rows_checked": rows_checked,
            "last_bad_prime": last_bad,
            "tail_start_prime": tail_start,
            "terminal_block_min_margin_lo": (
                mpmath.nstr(post_bad_min, dps) if post_bad_min is not None else None
            ),
            "terminal_block_min_margin_prime": post_bad_min_prime,
        }


def _entry_parameter(kappa_text: str, q: int, dps: int) -> tuple[dict[str, Any], Any]:
    """Certify eta_K(q), returning its JSON record and interval value."""
    iv = mpmath.iv
    with working_precision(dps), using_field(MpiField()):
        q_iv = mpmath.mpi(q)
        log_q = iv.log(q_iv)
        values: list[Any] = []
        best_index = 0
        # Double the search range until the saturation bound excludes the tail.
        candidate_limit = 2
        while True:
            values = [
                iv.log(rho_pp(q, a, kappa_text)) / (a * log_q)
                for a in range(1, candidate_limit + 1)
            ]
            point_values = [_point(v) for v in values]
            best_index = max(range(candidate_limit), key=point_values.__getitem__)
            best = values[best_index]
            tail_upper = iv.log(q_iv / (q_iv - 1)) / ((candidate_limit + 1) * log_q)
            if tail_upper.b < best.a:
                break
            candidate_limit *= 2
            if candidate_limit > 4096:
                raise RuntimeError("failed to bound the entry-parameter exponent")

        competitors_separated = all(
            values[i].b < best.a for i in range(candidate_limit) if i != best_index
        )
        if not competitors_separated:
            raise RuntimeError("entry-parameter maximizer not separated by intervals")

        return (
            {
                "prime": q,
                "exponent": best_index + 1,
                "epsilon_lo": _lo(best, dps),
                "epsilon_hi": _hi(best, dps),
                "candidate_limit": candidate_limit,
                "tail_eta_upper": _hi(tail_upper, dps),
            },
            best,
        )


def _local_cap(p: int, epsilon_lo: mpmath.mpf, kappa_text: str, dps: int) -> tuple[int, str]:
    """Find a certified exponent after which the local objective decreases."""
    iv = mpmath.iv
    with working_precision(dps), using_field(MpiField()):
        log_p = iv.log(p)
        for a in range(1, 4097):
            ratio_bound = sigma_minus_one_pp(p, a + 1) / sigma_minus_one_pp(p, a - 1)
            slope_bound = iv.log(ratio_bound) / log_p
            if slope_bound.b < epsilon_lo:
                return a, _hi(slope_bound, dps)
    raise RuntimeError(f"failed to find a local exponent cap for p={p}")


def _propose_hull(point_y: list[mpmath.mpf]) -> list[int]:
    """Propose upper-hull vertices for later directed-interval certification."""
    hull: list[int] = []
    for a, y in enumerate(point_y):
        while len(hull) >= 2:
            a0, a1 = hull[-2], hull[-1]
            s01 = (point_y[a1] - point_y[a0]) / (a1 - a0)
            s12 = (y - point_y[a1]) / (a - a1)
            if s01 <= s12:
                hull.pop()
            else:
                break
        hull.append(a)
    return hull


def _local_profile(
    p: int,
    kappa_text: str,
    epsilon_lo: mpmath.mpf,
    epsilon_hi: mpmath.mpf,
    boundary_prime: int,
    boundary_exponent: int,
    dps: int,
) -> tuple[dict[str, Any], list[_CertifiedEvent]]:
    """Certify one local upper hull and return its relevant transitions."""
    iv = mpmath.iv
    cap, tail_slope_upper = _local_cap(p, epsilon_lo, kappa_text, dps)

    with working_precision(dps * 2), using_field(MpfField()):
        log_p_point = mpmath.log(p)
        y_point = [mpmath.mpf(0)] + [
            mpmath.log(rho_pp(p, a, kappa_text)) / log_p_point for a in range(1, cap + 1)
        ]
        hull = _propose_hull(y_point)

    with working_precision(dps), using_field(MpiField()):
        log_p = iv.log(p)
        y = [mpmath.mpi(0)] + [iv.log(rho_pp(p, a, kappa_text)) / log_p for a in range(1, cap + 1)]
        event_intervals: list[Any] = []
        for left, right in zip(hull, hull[1:]):
            event_intervals.append((y[right] - y[left]) / (right - left))

        for first, second in zip(event_intervals, event_intervals[1:]):
            if first.a <= second.b:
                raise RuntimeError(f"local hull slopes not separated for p={p}")

        for left, right in zip(hull, hull[1:]):
            for a in range(left + 1, right):
                chord = y[left] + (a - left) * (y[right] - y[left]) / (right - left)
                if y[a].b >= chord.a:
                    raise RuntimeError(f"non-hull point not separated below chord for p={p}, a={a}")

        relevant: list[_CertifiedEvent] = []
        for (left, right), event in zip(zip(hull, hull[1:]), event_intervals):
            if event.a > epsilon_hi:
                relevant.append(_CertifiedEvent(p, left, right, event))
            elif event.b < epsilon_lo:
                continue
            elif p == boundary_prime and left == 0 and right == boundary_exponent:
                relevant.append(_CertifiedEvent(p, left, right, event, boundary_event=True))
            else:
                raise RuntimeError(
                    f"event at p={p}, {left}->{right} is not separated from boundary epsilon"
                )

        profile = {
            "prime": p,
            "exponent_cap": cap,
            "tail_slope_upper": tail_slope_upper,
            "hull_vertices": hull,
            "relevant_event_count": len(relevant),
        }
        return profile, relevant


def _endpoint_margin(
    factorization: dict[int, int],
    kappa_text: str,
    dps: int,
    log_prime_cache: dict[int, Any],
    log_rho_cache: dict[tuple[int, int], Any],
) -> tuple[Any, Any, Any]:
    """Return the endpoint margin, reusing previously computed local logarithms."""
    iv = mpmath.iv
    with working_precision(dps), using_field(MpiField()):
        log_n = mpmath.mpi(0)
        log_q = mpmath.mpi(0)
        for p, a in factorization.items():
            if a <= 0:
                continue
            log_p = log_prime_cache.get(p)
            if log_p is None:
                log_p = iv.log(p)
                log_prime_cache[p] = log_p
            log_n += a * log_p

            key = (p, a)
            log_rho = log_rho_cache.get(key)
            if log_rho is None:
                log_rho = iv.log(rho_pp(p, a, kappa_text))
                log_rho_cache[key] = log_rho
            log_q += log_rho
        log_robin_scale = iv.euler + iv.log(iv.log(log_n))
        return log_robin_scale - log_q, log_n, log_q


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


def _initial_integer_scan(kappa_text: str, upper_exclusive: int, dps: int) -> dict[str, Any]:
    """Certify every integer from 5041 through ``upper_exclusive - 1``."""
    iv = mpmath.iv
    if upper_exclusive <= 5041:
        return {
            "lower": 5041,
            "upper": 5040,
            "count": 0,
            "all_positive": True,
            "minimum_margin_lo": None,
            "minimum_at_n": None,
        }
    minimum: mpmath.mpf | None = None
    minimum_n: int | None = None
    spf = _smallest_prime_factor_table(upper_exclusive - 1)
    with working_precision(dps), using_field(MpiField()):
        e_gamma = iv.exp(iv.euler)
        kappa = mpmath.mpi(kappa_text)
        one = mpmath.mpi(1)
        root = one / (kappa - one)
        local_base_cache: dict[tuple[int, int], Any] = {}
        for n in range(5041, upper_exclusive):
            scale = e_gamma * iv.log(iv.log(n))
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
                    local_base = p_iv**a * rksigma_pp(p, a, kappa) / kappa
                    local_base_cache[key] = local_base
                base_product *= local_base
            margin = scale - base_product**root
            if margin.a <= 0:
                raise RuntimeError(f"initial integer scan failed at n={n}")
            lower = mpmath.mpf(margin.a)
            if minimum is None or lower < minimum:
                minimum = lower
                minimum_n = n
    return {
        "lower": 5041,
        "upper": upper_exclusive - 1,
        "count": upper_exclusive - 5041,
        "all_positive": True,
        "minimum_margin_lo": mpmath.nstr(minimum, dps) if minimum is not None else None,
        "minimum_at_n": minimum_n,
    }


def compute_certificate(kappa_text: str, dps: int = 60) -> dict[str, Any]:
    """Generate a bounded-parameter certificate at parameter ``kappa_text``."""
    with unbounded_decimal_integers(), working_precision(dps):
        return _compute_certificate_at_precision(kappa_text, dps)


def _compute_certificate_at_precision(kappa_text: str, dps: int) -> dict[str, Any]:
    """Generate the certificate assuming the requested precision is active."""
    kappa_point = mpmath.mpf(kappa_text)
    if not mpmath.isfinite(kappa_point) or kappa_point <= 1:
        raise ValueError("kappa must be a finite real number greater than 1")

    tail = _finite_tail_scan(kappa_text, dps)
    q = int(tail["tail_start_prime"])
    entry, entry_interval = _entry_parameter(kappa_text, q, dps)
    epsilon_lo = mpmath.mpf(entry_interval.a)
    epsilon_hi = mpmath.mpf(entry_interval.b)

    profiles: list[dict[str, Any]] = []
    events: list[_CertifiedEvent] = []
    for p in primerange(2, q + 1):
        profile, local_events = _local_profile(
            int(p),
            kappa_text,
            epsilon_lo,
            epsilon_hi,
            q,
            int(entry["exponent"]),
            dps,
        )
        profiles.append(profile)
        events.extend(local_events)

    # Point arithmetic proposes the hulls. Ordering uses the certified interval
    # transitions themselves, not their decimal renderings.
    events.sort(key=lambda event: mpmath.mpf(event.interval.a), reverse=True)
    for current, following in zip(events, events[1:]):
        if current.interval.a <= following.interval.b:
            raise RuntimeError(
                "global endpoint events are not separated: "
                f"{current.prime}:{current.from_exponent}->{current.to_exponent} and "
                f"{following.prime}:{following.from_exponent}->{following.to_exponent}"
            )

    factorization: dict[int, int] = {}
    endpoint_rows: list[dict[str, Any]] = []
    n = 1
    first_above_5040: int | None = None
    endpoint_min_margin: mpmath.mpf | None = None
    endpoint_min_n: int | None = None
    previous_n = 1
    endpoint_log_prime_cache: dict[int, Any] = {}
    endpoint_log_rho_cache: dict[tuple[int, int], Any] = {}

    for index, event in enumerate(events, start=1):
        p = event.prime
        old = event.from_exponent
        new = event.to_exponent
        if factorization.get(p, 0) != old:
            raise RuntimeError(f"event sequence state mismatch at p={p}: expected exponent {old}")
        n *= p ** (new - old)
        factorization[p] = new
        if n <= previous_n:
            raise RuntimeError("canonical endpoint sequence failed to increase")
        previous_n = n

        row: dict[str, Any] = {
            "index": index,
            "event_prime": p,
            "from_exponent": old,
            "to_exponent": new,
            "epsilon_lo": _lo(event.interval, dps),
            "epsilon_hi": _hi(event.interval, dps),
            "endpoint": str(n),
            "largest_prime": max(factorization),
        }
        if n > 5040:
            margin, log_n, _ = _endpoint_margin(
                factorization,
                kappa_text,
                dps,
                endpoint_log_prime_cache,
                endpoint_log_rho_cache,
            )
            if margin.a <= 0:
                raise RuntimeError(f"Robin endpoint margin failed at N={n}")
            lower = mpmath.mpf(margin.a)
            row["robin_log_margin_lo"] = _lo(margin, dps)
            row["robin_log_margin_hi"] = _hi(margin, dps)
            row["log_n_lo"] = _lo(log_n, dps)
            row["log_n_hi"] = _hi(log_n, dps)
            if first_above_5040 is None:
                first_above_5040 = n
            if endpoint_min_margin is None or lower < endpoint_min_margin:
                endpoint_min_margin = lower
                endpoint_min_n = n
        endpoint_rows.append(row)

    if not endpoint_rows:
        raise RuntimeError("no endpoint events were generated")
    final_event = events[-1]
    final = endpoint_rows[-1]
    if final_event.prime != q or final_event.from_exponent != 0:
        raise RuntimeError("final event is not the certified tail-prime entry event")
    if first_above_5040 is None:
        raise RuntimeError("endpoint sequence did not reach beyond 5040")

    initial_scan = _initial_integer_scan(kappa_text, first_above_5040, dps)

    return {
        "schema_version": "bounded-parameter-endpoints-1.0.0",
        "status": "directed-interval-finite-certificate",
        "precision_dps": dps,
        "kappa": str(kappa_text),
        "theorem_scope": f"all n > 5040 and all real 1 < kappa <= {kappa_text}",
        "finite_tail": tail,
        "entry_parameter": entry,
        "local_profiles": profiles,
        "global_events": [event.record(dps) for event in events],
        "endpoint_verification": {
            "event_count": len(events),
            "first_endpoint_above_5040": str(first_above_5040),
            "tail_endpoint": final["endpoint"],
            "all_positive": True,
            "minimum_log_margin_lo": (
                mpmath.nstr(endpoint_min_margin, dps) if endpoint_min_margin is not None else None
            ),
            "minimum_log_margin_endpoint": (
                str(endpoint_min_n) if endpoint_min_n is not None else None
            ),
            "rows": endpoint_rows,
        },
        "initial_integer_scan": initial_scan,
    }


def main(
    argv: Sequence[str] | None = None,
    *,
    default_output_dir: Path | None = None,
    default_kappa: str | None = None,
) -> int:
    """Generate one bounded-parameter certificate from the command line."""
    parser = argparse.ArgumentParser(
        description="Generate a directed-interval bounded-parameter certificate."
    )
    parser.add_argument(
        "--kappa",
        default=default_kappa,
        required=default_kappa is None,
        help=(
            "fixed endpoint parameter K > 1"
            + (f" (default: {default_kappa})" if default_kappa is not None else "")
        ),
    )
    parser.add_argument("--dps", type=int, default=60, help="interval precision in decimal digits")
    parser.add_argument("--output", type=Path, help="output JSON path")
    args = parser.parse_args(argv)
    output = args.output or default_output_path(args.kappa, default_output_dir)
    result = compute_certificate(args.kappa, args.dps)
    write_json(output, result)
    print(f"Wrote {output}")
    print(f"Tail starts at P={result['finite_tail']['tail_start_prime']}")
    print(f"Canonical events: {result['endpoint_verification']['event_count']}")
    print(
        "First endpoint above 5040: "
        f"{result['endpoint_verification']['first_endpoint_above_5040']}"
    )
    print(
        "Minimum endpoint log margin: "
        f"{result['endpoint_verification']['minimum_log_margin_lo']}"
    )
    print("Initial direct scan minimum at n=" f"{result['initial_integer_scan']['minimum_at_n']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
