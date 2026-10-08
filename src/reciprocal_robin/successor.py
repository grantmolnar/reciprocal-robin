"""Enumerate reciprocal kappa-CA candidates by threshold events.

As ``eps`` decreases, one prime exponent changes at a time. The event taking
prime ``p`` to exponent ``ell`` occurs at ``F^[k](p, ell)``. A max-heap keeps
the next event for each active prime and the first event for the next prime.

Point arithmetic gives the numerical event order used for exploration. With
interval arithmetic, heap priorities are interval midpoints and therefore do
not certify the order of near ties.
"""

from heapq import heappop, heappush
from typing import Iterator

from sympy import nextprime

from ._field import _to_float, current_field
from ._types import NumericInput, Scalar
from ._validation import require_kappa
from .arithmetic import F_local
from .factorization import Factorization


class CASuccessorEnumerator:
    """Yield reciprocal-CA events in decreasing numerical ``eps``.

    Each item is ``(eps, factorization)`` after one prime exponent changes. The
    numerical backend is fixed at construction because the heap stores values
    from that field.
    """

    def __init__(self, k: NumericInput) -> None:
        require_kappa(k)
        self.k = k
        self._field_name = current_field().name
        self._exponents: dict[int, int] = {}
        # heapq is a min-heap, so negated priorities give max-heap behavior.
        # (p, ell) provides deterministic tie-breaking after the priority.
        self._heap: list[tuple[float, int, int, Scalar]] = []
        self._push_event(2, 1)

    def _require_same_field(self) -> None:
        current_name = current_field().name
        if current_name != self._field_name:
            raise RuntimeError(
                "CASuccessorEnumerator cannot change numerical backend while "
                f"iterating (created under {self._field_name}, now {current_name})"
            )

    def _push_event(self, p: int, ell: int) -> None:
        self._require_same_field()
        eps = F_local(p, ell, self.k)
        # Interval events have no total order; midpoint priority is only a
        # numerical ordering device and cannot certify separation.
        priority = -_to_float(eps)
        heappush(self._heap, (priority, p, ell, eps))

    def current_factorization(self) -> Factorization:
        """Return the current candidate factorization as an immutable value."""
        return Factorization._from_trusted_dict(self._exponents)

    def __iter__(self) -> Iterator[tuple[Scalar, Factorization]]:
        return self

    def __next__(self) -> tuple[Scalar, Factorization]:
        self._require_same_field()
        if not self._heap:
            raise StopIteration
        _, p, ell, eps = heappop(self._heap)

        previous = self._exponents.get(p, 0)
        assert ell == previous + 1, (
            f"CA event ({p}, {ell}) out of order: prime {p} is at "
            f"exponent {previous}, cannot jump to {ell}"
        )
        self._exponents[p] = ell

        self._push_event(p, ell + 1)
        if previous == 0:
            self._push_event(int(nextprime(p)), 1)

        return eps, self.current_factorization()
