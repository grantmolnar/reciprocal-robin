"""Lazy reciprocal-CA walks over explicit epsilon grids.

Each :class:`CAStep` records one reciprocal-CA factorization, its normalized
Robin ratio, and the verdict supplied by the active numerical field.
:mod:`.ca_grid` provides presentation and serialization.
"""

from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from mpmath import mpf

from ._field import MpfField, Trilean, current_field, using_field
from ._types import NumericInput, Scalar
from ._validation import require_kappa
from .ca import reciprocal_ca
from .factorization import Factorization


@dataclass(frozen=True)
class CAStep:
    """One reciprocal-CA sample with structural data and a field verdict."""

    eps: Scalar
    factorization: Factorization
    ratio: Scalar
    verdict: Trilean

    @property
    def omega(self) -> int:
        return self.factorization.omega

    @property
    def p_max(self) -> int:
        return self.factorization.p_max

    @property
    def log_n(self) -> Scalar:
        return self.factorization.log_n()


def walk_ca(k: NumericInput, eps_values: Iterable[NumericInput]) -> Iterator[CAStep]:
    """Yield :class:`CAStep` records for the supplied epsilon values.

    Candidate factorizations are branch-dependent and are therefore constructed
    with point arithmetic. Once a factorization is fixed, its ratio is evaluated
    in the field active when iteration began. Thus an :class:`MpiField` caller
    obtains interval ratios without using intervals to choose the factorization.

    ``CAStep.eps`` is the ``mpf`` value used by the constructor. Empty
    factorizations are skipped.
    """
    require_kappa(k)
    evaluation_field = current_field()
    evaluation_k = evaluation_field.number(k)
    point_k = mpf(k)

    for eps in eps_values:
        point_eps = mpf(eps)
        with using_field(MpfField()):
            fac = reciprocal_ca(point_k, point_eps)
        if not fac:
            continue
        # Factorization methods read current_field(); restore the field captured
        # when the iterator was created before evaluating each fixed candidate.
        with using_field(evaluation_field):
            ratio = fac.robin_ratio(evaluation_k)
            verdict = evaluation_field.is_strictly_less_than_one(ratio)
        yield CAStep(eps=point_eps, factorization=fac, ratio=ratio, verdict=verdict)
