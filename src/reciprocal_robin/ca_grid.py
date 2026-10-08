"""Evaluate reciprocal-CA numbers on explicit epsilon grids.

Point scans report sampled ratios. Interval audits classify the normalized
Robin ratio rigorously at each sampled reciprocal-CA number. Neither operation
controls integers between the sampled candidates.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import asdict, dataclass
from datetime import datetime, timezone

import mpmath

from ._field import MpiField, Trilean, _to_float, current_field, using_field
from ._precision import working_precision
from ._serialization import json_text
from ._types import NumericInput, Scalar
from ._validation import require_integer, require_kappa
from .factorization import Factorization
from .walker import CAStep, walk_ca

StepObserver = Callable[[Scalar, Factorization, Scalar], None]


def max_robin_ratio_on_ca_grid(
    k: NumericInput,
    eps_values: Iterable[NumericInput],
    observer: StepObserver | None = None,
) -> Scalar:
    """Return the largest point Robin ratio on an explicit reciprocal-CA grid.

    No interpolation or tail estimate is inferred from the grid. Use
    :func:`audit_ca_grid` for directed-interval verdicts at the sampled points.
    """
    require_kappa(k)
    field = current_field()
    if field.name != "mpf":
        raise TypeError("max_robin_ratio_on_ca_grid is point-valued; use MpfField explicitly")
    eps_inputs = tuple(eps_values)
    if not eps_inputs:
        raise ValueError("eps_values must contain at least one value")

    max_ratio: Scalar | None = None
    for step in walk_ca(k, eps_inputs):
        if max_ratio is None or step.ratio > max_ratio:
            max_ratio = step.ratio
        if observer is not None:
            observer(step.eps, step.factorization, step.ratio)
    if max_ratio is None:
        raise ValueError("epsilon grid produced no nontrivial reciprocal-CA candidates")
    return max_ratio


def print_step(eps: Scalar, fac: Factorization, ratio: Scalar) -> None:
    """Print one row suitable as an observer for a point grid scan."""
    print(
        f"eps={float(eps):+.4f}  omega={fac.omega:6d}  P*={fac.p_max:10d}  "
        f"log n={_to_float(fac.log_n()):12.2f}  ratio={_to_float(ratio):.6f}"
    )


@dataclass(frozen=True)
class CertifiedStep:
    """Serializable directed-interval verdict for one sampled CA point."""

    eps: str
    omega: int
    p_max: int
    log_n_lo: str
    log_n_hi: str
    ratio_lo: str
    ratio_hi: str
    verdict: str


@dataclass(frozen=True)
class CAWalkReport:
    """Directed-interval verdicts for reciprocal-CA points on an explicit grid.

    The report concerns only the listed points; it says nothing about
    intervening integers, unsampled epsilon values, or an unbounded tail.
    ``steps`` is immutable.
    """

    kappa: str
    precision_dps: int
    steps: tuple[CertifiedStep, ...]
    all_steps_resolved: bool
    ambiguity_witness: CertifiedStep | None
    timestamp: str
    schema_version: str = "ca-walk-report-1.1.0"

    def to_json(self) -> str:
        """Serialize the report as human-readable JSON."""
        return json_text(asdict(self))


def _interval_endpoints(value: Scalar) -> tuple[str, str]:
    """Serialize the lower and upper endpoints of a point or interval value."""
    if hasattr(value, "a") and hasattr(value, "b"):
        lo = mpmath.nstr(mpmath.mpf(value.a), mpmath.mp.dps)
        hi = mpmath.nstr(mpmath.mpf(value.b), mpmath.mp.dps)
    else:
        lo = hi = mpmath.nstr(value, mpmath.mp.dps)
    return lo, hi


def _serialize_step(step: CAStep) -> CertifiedStep:
    ratio_lo, ratio_hi = _interval_endpoints(step.ratio)
    log_n_lo, log_n_hi = _interval_endpoints(step.log_n)
    return CertifiedStep(
        eps=mpmath.nstr(step.eps, mpmath.mp.dps),
        omega=step.omega,
        p_max=step.p_max,
        log_n_lo=log_n_lo,
        log_n_hi=log_n_hi,
        ratio_lo=ratio_lo,
        ratio_hi=ratio_hi,
        verdict=step.verdict.name,
    )


def _resolve_verdict_with_bump(
    ratio: Scalar,
    fac: Factorization,
    k: NumericInput,
    bumped_dps: int,
) -> tuple[Scalar, Trilean]:
    """Classify one ratio, with one higher-precision retry if it straddles one."""
    field = current_field()
    verdict = field.is_strictly_less_than_one(ratio)
    if verdict is Trilean.AMBIGUOUS:
        with working_precision(bumped_dps):
            # Reparse the original input so textual kappa values gain the bumped
            # precision instead of inheriting a lower-precision mpf rounding.
            bumped_k = field.number(k)
            ratio = fac.robin_ratio(bumped_k)
            verdict = field.is_strictly_less_than_one(ratio)
    return ratio, verdict


def audit_ca_grid(
    k: NumericInput,
    eps_values: Iterable[NumericInput],
    *,
    base_dps: int = 50,
    bumped_dps: int = 200,
) -> CAWalkReport:
    """Classify the Robin ratio at each CA point on an explicit epsilon grid.

    Text inputs are parsed at ``base_dps``. If a ratio interval straddles one,
    the fixed factorization is evaluated once more at ``bumped_dps``. An
    already-rounded ``mpf`` input cannot recover digits when the precision is
    raised.
    """
    require_integer(base_dps, "base_dps", minimum=2)
    require_integer(bumped_dps, "bumped_dps", minimum=2)
    if bumped_dps < base_dps:
        raise ValueError("bumped_dps must be at least base_dps")

    eps_inputs = tuple(eps_values)
    if not eps_inputs:
        raise ValueError("eps_values must contain at least one value")

    steps: list[CertifiedStep] = []
    ambiguity_witness: CertifiedStep | None = None
    all_steps_resolved = True

    with working_precision(base_dps):
        require_kappa(k)
        kappa_text = mpmath.nstr(mpmath.mpf(k), mpmath.mp.dps)
        with using_field(MpiField()):
            for step in walk_ca(k, eps_inputs):
                ratio, verdict = _resolve_verdict_with_bump(
                    step.ratio, step.factorization, k, bumped_dps
                )
                resolved = CAStep(
                    eps=step.eps,
                    factorization=step.factorization,
                    ratio=ratio,
                    verdict=verdict,
                )
                serialized = _serialize_step(resolved)
                steps.append(serialized)
                if verdict is Trilean.AMBIGUOUS:
                    all_steps_resolved = False
                    ambiguity_witness = serialized
                    break

    if not steps:
        raise ValueError("epsilon grid produced no nontrivial reciprocal-CA candidates")

    return CAWalkReport(
        kappa=kappa_text,
        precision_dps=base_dps,
        steps=tuple(steps),
        all_steps_resolved=all_steps_resolved,
        ambiguity_witness=ambiguity_witness,
        timestamp=datetime.now(timezone.utc).isoformat(),
    )
