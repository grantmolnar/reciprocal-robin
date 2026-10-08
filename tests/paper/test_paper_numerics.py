"""Regression tests for manuscript numerics and registry integrity.

The numerical claims printed in "Reciprocal LCM powers and Robin's criterion"
(Molnar 2026) are registered in :mod:`reciprocal_robin.paper_numerics`. This
module checks the complete registry and retains focused structural regressions
for Example 5.6. A failure means either the implementation drifted or the
manuscript and its numerical registry need to be reconciled.

The direct Example 5.6 regressions pin ``MpfField`` because they compare
computed values with printed paper numbers via ``float(...)``. The MpiField
default is exercised separately in ``test_field.py`` and through the
interval-arithmetic walker paths in ``test_ca_grid.py``.
"""

from mpmath import mpf
import pytest

from reciprocal_robin import MpfField, using_field
from reciprocal_robin.exploration import reciprocal_ca
from reciprocal_robin.paper_numerics import CLAIMS, Status, verify_all


@pytest.fixture(autouse=True)
def _pin_mpf_field():
    """Pin MpfField for all tests in this module (paper-number baselines)."""
    with using_field(MpfField()):
        yield


# ---- Example 5.6: structural columns at k = 2, eps = -0.40 ------------------


def test_example_5_6_structure():
    """omega = 11 and P_max = 31 at the (k=2, eps=-0.40) reciprocal CA number.

    These are the structural quantities that *do* match the printed table.
    """
    fac = reciprocal_ca(mpf(2), mpf("-0.40"))
    assert fac.omega == 11
    assert fac.p_max == 31


def test_example_5_6_factorization_and_ratio():
    """The threshold rule gives exponent 2, not 3, at the prime 2."""
    fac = reciprocal_ca(mpf(2), mpf("-0.40"))
    assert fac.primes_to_exponents[2] == 2
    assert all(exponent == 1 for prime, exponent in fac.primes_to_exponents.items() if prime != 2)
    ratio = fac.robin_ratio(mpf(2))
    assert abs(float(ratio) - 0.46256784429434305) < 1e-14


# ---- Registered manuscript numerics ------------------------------------------


def test_every_registered_claim_is_confirmed():
    verdicts = verify_all()
    failures = [v for v in verdicts if v.status is not Status.CONFIRMED]
    assert not failures, [(v.claim.locator, v.status, v.note) for v in failures]


def test_registry_is_immutable_and_every_claim_has_distinct_locator():
    assert isinstance(CLAIMS, tuple)
    locators = [claim.locator for claim in CLAIMS]
    assert len(locators) == len(set(locators))


def test_withdrawn_transfer_constants_are_not_registered():
    text = " ".join(claim.description + " " + claim.locator for claim in CLAIMS)
    assert "42055" not in text
    assert "kappa_MP" not in text
