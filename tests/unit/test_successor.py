"""Successor-enumerator tests and the CA-endpoint counterexample."""

import math
from itertools import islice

import pytest
from mpmath import mp, mpf

from reciprocal_robin import robin_ratio
from reciprocal_robin.successor import CASuccessorEnumerator


@pytest.fixture(autouse=True)
def _pin_mpf_field():
    from reciprocal_robin import MpfField, using_field

    with using_field(MpfField()):
        yield


def test_successor_first_event_matches_F_local():
    from reciprocal_robin import F_local

    enum = CASuccessorEnumerator(mpf(2))
    first_eps, first_fac = next(enum)
    expected = F_local(2, 1, mpf(2))
    assert mp.almosteq(first_eps, expected, rel_eps=mpf("1e-30"), abs_eps=mpf("1e-30"))
    assert first_fac.primes_to_exponents == {2: 1}


def test_ca_endpoint_local_maximum_claim_is_false():
    """CA endpoints alone do not control every intervening integer."""
    enum = CASuccessorEnumerator(mpf(2))
    first_four = list(islice(enum, 4))
    values = [math.prod(p**a for p, a in fac.primes_to_exponents.items()) for _, fac in first_four]
    assert values == [2, 6, 30, 210]
    r30 = float(robin_ratio(30, mpf(2)))
    r32 = float(robin_ratio(32, mpf(2)))
    r210 = float(robin_ratio(210, mpf(2)))
    assert r32 > max(r30, r210)


def test_successor_rejects_backend_change_midstream():
    from reciprocal_robin import MpiField, using_field

    enum = CASuccessorEnumerator(mpf(2))
    next(enum)
    with using_field(MpiField()):
        with pytest.raises(RuntimeError, match="cannot change numerical backend"):
            next(enum)
