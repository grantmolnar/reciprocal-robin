"""Update the normalized Robin ratio along a reciprocal-CA successor walk.

Consecutive candidates differ at one prime exponent. The local factor at that
prime and ``log N`` can therefore be updated without recomputing the full
product. The walker accepts any iterable of ``(epsilon, factorization)`` pairs;
the update identity does not depend on a particular event generator.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from ._field import current_field
from ._types import NumericInput, Scalar
from ._validation import require_kappa
from .arithmetic import rksigma_pp
from .factorization import Factorization

IncrementalStep = tuple[Scalar, Factorization, Scalar]
SuccessorStep = tuple[Scalar, Factorization]


class SuccessorInvariantError(RuntimeError):
    """Raised when consecutive factorizations violate the one-event contract."""


@dataclass(frozen=True)
class IncrementalRatioState:
    """Read-only diagnostic snapshot of the running numerical state."""

    product: Scalar
    log_n: Scalar
    exponents: tuple[tuple[int, int], ...]


class IncrementalRatioWalker:
    """Track the Robin ratio along a one-event successor stream.

    Output begins at the first endpoint where the real Robin ratio is defined;
    the structural ``n=2`` event is consumed internally. The numerical backend
    is fixed at construction because the running product and logarithm are field
    elements.
    """

    def __init__(self, successors: Iterable[SuccessorStep], k: NumericInput) -> None:
        require_kappa(k)
        self._successors = iter(successors)
        self._k = k
        field = current_field()
        self._field_name = field.name
        self._product = field.number(1)
        self._log_n = field.number(0)
        self._exponents_seen: dict[int, int] = {}

    def _require_same_field(self) -> None:
        current_name = current_field().name
        if current_name != self._field_name:
            raise RuntimeError(
                "IncrementalRatioWalker cannot change numerical backend while "
                f"iterating (created under {self._field_name}, now {current_name})"
            )

    @property
    def state(self) -> IncrementalRatioState:
        """Return the current product, logarithm, and exponent vector."""
        return IncrementalRatioState(
            product=self._product,
            log_n=self._log_n,
            exponents=tuple(sorted(self._exponents_seen.items())),
        )

    def __iter__(self) -> Iterator[IncrementalStep]:
        return self

    def __next__(self) -> IncrementalStep:
        """Advance to the next endpoint on which the real Robin ratio is defined."""
        self._require_same_field()
        field = current_field()
        k = self._k

        while True:
            eps, factorization = next(self._successors)
            new_pairs = factorization.primes_to_exponents
            self._check_existing_primes(new_pairs)
            event_prime, event_exponent = self._identify_single_event(new_pairs)

            self._log_n = self._log_n + field.log(field.number(event_prime))

            previous_exponent = event_exponent - 1
            p_field = field.number(event_prime)
            new_local = p_field**event_exponent * rksigma_pp(event_prime, event_exponent, k) / k
            if previous_exponent == 0:
                old_local = field.number(1)
            else:
                old_local = (
                    p_field**previous_exponent * rksigma_pp(event_prime, previous_exponent, k) / k
                )
            self._product = self._product * new_local / old_local
            self._exponents_seen[event_prime] = event_exponent

            if self._log_n <= field.number(1):
                continue

            rho = self._product ** (field.number(1) / (k - 1))
            ratio = rho / (field.exp(field.euler_gamma()) * field.log(self._log_n))
            return eps, factorization, ratio

    def _check_existing_primes(self, new_pairs: dict[int, int]) -> None:
        """Check that previously present prime exponents do not decrease."""
        for prime, previous in self._exponents_seen.items():
            current = new_pairs.get(prime)
            if current is None:
                raise SuccessorInvariantError(f"Prime {prime} disappeared between CA steps")
            if current < previous:
                raise SuccessorInvariantError(
                    f"Prime {prime} exponent decreased: {previous} -> {current}"
                )

    def _identify_single_event(self, new_pairs: dict[int, int]) -> tuple[int, int]:
        """Return the unique ``(prime, exponent)`` increment in a successor step."""
        event: tuple[int, int] | None = None
        for prime, exponent in new_pairs.items():
            previous = self._exponents_seen.get(prime, 0)
            if exponent == previous:
                continue
            if event is not None:
                raise SuccessorInvariantError(
                    "More than one prime changed between CA steps: " f"{event[0]} and {prime}"
                )
            if exponent != previous + 1:
                raise SuccessorInvariantError(
                    f"Prime {prime} advanced by more than one exponent: "
                    f"{previous} -> {exponent}"
                )
            event = (prime, exponent)

        if event is None:
            raise SuccessorInvariantError(
                "No prime changed between CA steps; successor stream violated "
                "the one-event-per-step invariant"
            )
        return event
