"""Registry for numerical values printed in the manuscript.

Each :class:`Claim` records the printed value, its comparison tolerance, and a
producer that recomputes it from the numerical implementation. The registry is
immutable after import. The theorem certificates under
:mod:`reciprocal_robin.certification` are verified separately.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum

import mpmath
from mpmath import mpf

from ._field import MpfField, using_field
from ._precision import working_precision
from ._types import Scalar
from ._validation import require_precision

Producer = Callable[[], Scalar]


class Status(Enum):
    """Outcome of comparing one recomputed value with its registered value."""

    CONFIRMED = "confirmed"
    SLIP = "slip"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class Claim:
    """One manuscript numerical value and its reproducibility contract."""

    locator: str
    description: str
    expected: str
    tolerance: str
    produce: Producer


@dataclass(frozen=True)
class Verdict:
    """Result of recomputing one :class:`Claim`."""

    claim: Claim
    computed: Scalar | None
    status: Status
    note: str = ""


def verify(claim: Claim, dps: int = 80) -> Verdict:
    """Recompute one claim at ``dps`` decimal digits and compare its tolerance."""
    require_precision(dps)
    with working_precision(dps):
        try:
            with using_field(MpfField()):
                value = claim.produce()
        except Exception as exc:  # one failed producer should not suppress the remaining claims
            return Verdict(claim, None, Status.UNKNOWN, f"exception: {exc}")

        if value is None:
            return Verdict(claim, None, Status.UNKNOWN, "producer returned None")

        value = mpf(value)
        delta = abs(value - mpf(claim.expected))
        if delta <= mpf(claim.tolerance):
            return Verdict(claim, value, Status.CONFIRMED)
        return Verdict(
            claim,
            value,
            Status.SLIP,
            f"computed {mpmath.nstr(value, 12)} vs expected {claim.expected} "
            f"(|delta|={mpmath.nstr(delta, 5)})",
        )


def verify_all(dps: int = 80) -> list[Verdict]:
    """Verify every registered manuscript value in registry order."""
    return [verify(claim, dps=dps) for claim in CLAIMS]


def _F_local(p: int, a: int) -> Producer:
    from .arithmetic import F_local

    return lambda: F_local(p, a, mpf(2))


def _ratio(n: int, k: str) -> Producer:
    from .normalized import robin_ratio

    return lambda: robin_ratio(n, mpf(k))


def _kappa_five_epsilon0() -> Producer:
    from .normalized_ca import rho_pp

    return lambda: mpmath.log(rho_pp(3, 3, 5) / rho_pp(3, 2, 5)) / mpmath.log(3)


def _kappa_five_low_support_margin() -> Producer:
    return lambda: (mpmath.exp(mpmath.euler) * mpmath.log(20 * mpmath.log(2)) - mpf(35) / 8)


def _kappa_five_p11_margin() -> Producer:
    return lambda: (
        mpmath.exp(mpmath.euler) * mpmath.log(2 * mpmath.log(2 * 3 * 5 * 7 * 11)) - mpf(77) / 16
    )


def _theta_19421_margin() -> Producer:
    from sympy import primerange

    def produce() -> Scalar:
        theta = sum((mpmath.log(p) for p in primerange(2, 19_422)), mpf(0))
        return mpf(19_183) - theta

    return produce


def _mertens_required_correction_p24137() -> Producer:
    from sympy import primerange

    def produce() -> Scalar:
        p = 24_137
        log_product = mpf(0)
        for q in primerange(2, p + 1):
            log_product -= mpmath.log1p(-mpf(1) / q)
        log_p = mpmath.log(p)
        return log_p**2 * (mpmath.exp(log_product - mpmath.euler) / log_p - 1)

    return produce


def _mertens_transfer_coefficient() -> Producer:
    return lambda: mpmath.exp(mpmath.euler) * (mpf("0.10952") + 1 / mpf("7.92"))


def _schoenfeld_log_comparison_margin() -> Producer:
    def produce() -> Scalar:
        log_p = mpmath.log(19_421)
        return mpmath.log(1 - 1 / (8 * log_p)) + 1 / (mpf("7.92") * log_p)

    return produce


def _bounded_tail_exp_factor() -> Producer:
    return lambda: mpmath.exp(mpf("0.010952"))


def _ca_metric(eps: str, metric: str) -> Producer:
    from .ca import reciprocal_ca

    def produce() -> Scalar:
        factorization = reciprocal_ca(mpf(2), mpf(eps))
        if metric == "omega":
            return mpf(factorization.omega)
        if metric == "p_max":
            return mpf(factorization.p_max)
        if metric == "ratio":
            return factorization.robin_ratio(mpf(2))
        raise ValueError(f"unknown CA metric: {metric}")

    return produce


def _build_claims() -> tuple[Claim, ...]:
    """Construct the fixed manuscript registry."""
    claims = [
        Claim(
            "CA table F(2,1)",
            "F^[2](2,1)",
            "0.16096404744368117",
            "1e-15",
            _F_local(2, 1),
        ),
        Claim(
            "CA table F(2,2)",
            "F^[2](2,2)",
            "-0.310744188373135",
            "1e-14",
            _F_local(2, 2),
        ),
    ]

    for eps, omega, p_max, ratio in (
        ("-0.30", 3, 5, "0.7357685847879101"),
        ("-0.35", 5, 11, "0.6149460620072547"),
        ("-0.40", 11, 31, "0.46256784429434305"),
        ("-0.43", 34, 139, "0.4854269606205743"),
        ("-0.45", 172, 1021, "0.4020402469891397"),
    ):
        claims.extend(
            [
                Claim(
                    f"CA example eps={eps} omega",
                    "reciprocal-CA number of prime factors",
                    str(omega),
                    "0",
                    _ca_metric(eps, "omega"),
                ),
                Claim(
                    f"CA example eps={eps} p_max",
                    "reciprocal-CA largest prime",
                    str(p_max),
                    "0",
                    _ca_metric(eps, "p_max"),
                ),
                Claim(
                    f"CA example eps={eps} ratio",
                    "reciprocal-CA Robin ratio",
                    ratio,
                    "5e-15",
                    _ca_metric(eps, "ratio"),
                ),
            ]
        )

    claims.extend(
        [
            Claim(
                "ratio(30,2)",
                "CA counterexample left endpoint",
                "0.7357685847879101",
                "1e-15",
                _ratio(30, "2"),
            ),
            Claim(
                "ratio(32,2)",
                "CA counterexample interior value",
                "0.8822741203451193",
                "1e-15",
                _ratio(32, "2"),
            ),
            Claim(
                "ratio(210,2)",
                "CA counterexample right endpoint",
                "0.5755898532050233",
                "1e-15",
                _ratio(210, "2"),
            ),
            Claim(
                "normalized-CA epsilon0",
                "supporting parameter for the endpoint 10800 at kappa=5",
                "0.04664156132118013343244854552401679038",
                "1e-35",
                _kappa_five_epsilon0(),
            ),
            Claim(
                "normalized-CA ratio(10800,5)",
                "Robin ratio at the first usable normalized endpoint",
                "0.8574354402857558",
                "1e-15",
                _ratio(10800, "5"),
            ),
            Claim(
                "kappa-five P=11 margin",
                "minimum finite squared-primorial Mertens margin",
                "0.06798576363532519141145628739529801948",
                "1e-35",
                _kappa_five_p11_margin(),
            ),
            Claim(
                "kappa-five low-support margin",
                "exp(gamma) log(20 log 2) - 35/8",
                "0.3078300704250456959372041778770399112",
                "1e-35",
                _kappa_five_low_support_margin(),
            ),
            Claim(
                "primorial threshold margin",
                "19183 - theta(19421)",
                "0.6540204126143589550689802752315266950",
                "1e-37",
                _theta_19421_margin(),
            ),
            Claim(
                "Mertens correction at P=24137",
                "required correction coefficient at the worst certified prime",
                "0.1095103735062916977688707829985333535",
                "1e-39",
                _mertens_required_correction_p24137(),
            ),
            Claim(
                "Mertens transfer coefficient",
                "exp(gamma) (0.10952 + 1/7.92)",
                "0.4199459322776549158224585080474472372",
                "2e-39",
                _mertens_transfer_coefficient(),
            ),
            Claim(
                "Schoenfeld logarithmic margin",
                "log(1 - 1/(8L)) + 1/(7.92L) at L=log(19421)",
                "0.00004705985536922150763093423104481860359",
                "1e-41",
                _schoenfeld_log_comparison_margin(),
            ),
            Claim(
                "bounded-parameter tail exponential factor",
                "exp(0.010952)",
                "1.011012192694765526109899904835059585969719802862363477613214",
                "1e-60",
                _bounded_tail_exp_factor(),
            ),
        ]
    )
    return tuple(claims)


CLAIMS: tuple[Claim, ...] = _build_claims()
