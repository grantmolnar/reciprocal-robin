"""Property tests for structural statements used in the current paper.

They check pointwise consequences of the power-mean representation and the
finite exceptional-set definition.
"""

from __future__ import annotations

import pytest
from hypothesis import given, settings, strategies as st
from mpmath import mpf
from sympy import divisor_sigma

from reciprocal_robin import MpfField, exceptional_set, rho_k, using_field
from reciprocal_robin.arithmetic import r_kappa_pa


@pytest.fixture(autouse=True)
def _pin_mpf_field():
    with using_field(MpfField()):
        yield


@settings(max_examples=100, deadline=None)
@given(
    n=st.integers(min_value=2, max_value=500),
    k_lo=st.floats(min_value=1.01, max_value=20.0, allow_nan=False, allow_infinity=False),
    delta=st.floats(min_value=0.01, max_value=20.0, allow_nan=False, allow_infinity=False),
)
def test_rho_k_is_strictly_increasing_and_below_sigma_minus_one(n, k_lo, delta):
    """For fixed n > 1, rho^[k](n) increases to sigma_{-1}(n)."""
    k1 = mpf(k_lo)
    k2 = mpf(k_lo + delta)
    rho1 = rho_k(n, k1)
    rho2 = rho_k(n, k2)
    sigma_minus_one = mpf(int(divisor_sigma(n))) / n
    assert rho1 < rho2 < sigma_minus_one


@settings(max_examples=100, deadline=None)
@given(
    p=st.sampled_from([2, 3, 5, 7, 11, 13]),
    a=st.integers(min_value=1, max_value=5),
    k=st.floats(min_value=1.01, max_value=50.0, allow_nan=False, allow_infinity=False),
)
def test_local_gap_factor_increases_with_prime_power_exponent(p, a, k):
    """The local gap factor increases with the prime-power exponent."""
    kk = mpf(k)
    assert r_kappa_pa(p, a, kk) < r_kappa_pa(p, a + 1, kk) < 1


@settings(max_examples=40, deadline=None)
@given(
    k_lo=st.sampled_from([mpf("1.5"), mpf(2), mpf(3), mpf(10)]),
    delta=st.sampled_from([mpf("0.5"), mpf(1), mpf(5)]),
    n_max=st.integers(min_value=10, max_value=80),
)
def test_finite_exceptional_set_grows_with_kappa(k_lo, delta, n_max):
    """On a fixed finite range, increasing kappa can only add exceptions."""
    low = set(exceptional_set(k_lo, n_max=n_max))
    high = set(exceptional_set(k_lo + delta, n_max=n_max))
    assert low.issubset(high)
