"""Property-based tests using Hypothesis.

These tests check structural identities and canonical representations in the
reciprocal sigma family. They are stronger than fixed numerical regressions
because they sample many small inputs.

Tests in this module pin ``MpfField`` because structural equalities
(``a == b`` on Dirichlet convolutions, e.g.) require exact comparison,
which interval arithmetic does not supply. The corresponding
properties under ``MpiField`` amount to interval containment, which is
covered by the walker/certify test paths.
"""

import pytest
from hypothesis import HealthCheck, given, settings, strategies as st
from mpmath import mp, mpf
from sympy import divisors, mobius, primerange

from reciprocal_robin import (
    Factorization,
    MpfField,
    rksigma_pp,
    sigma_minus_one_pp,
    using_field,
)
from reciprocal_robin.normalized import rho_k as rho_k_pointwise


@pytest.fixture(autouse=True)
def _pin_mpf_field():
    """Pin MpfField for property-based structural equality tests."""
    with using_field(MpfField()):
        yield


# Strategies -----------------------------------------------------------------

small_primes = st.sampled_from(list(primerange(2, 30)))
small_exponents = st.integers(min_value=1, max_value=4)
nice_kappas = st.sampled_from([2, 3, 4, 5, mpf("1.5"), mpf("2.5")])

slow = settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])


# Property 1: rksigma_pp matches its Dirichlet-convolution definition --------


@slow
@given(p=small_primes, a=small_exponents, k=nice_kappas)
def test_rksigma_pp_matches_dirichlet_convolution(p, a, k):
    """rksigma_pp(p, a, k) should equal sum_{d | p^a} mu(p^a/d) sigma_{-1}(d)^k.

    This is the defining identity of mu * sigma_{-1}^k as a Dirichlet
    convolution, computed two completely different ways.
    """
    k_mpf = mpf(k)
    formula_value = rksigma_pp(p, a, k_mpf)
    convolution_value = mpf(0)
    for d in divisors(p**a):
        mu = mobius(p**a // d)
        if mu == 0:
            continue
        # sigma_{-1}(d) = sigma(d) / d, computed via the prime-power formula
        # because d itself is a power of p (possibly d = 1).
        # Just use sigma_minus_one_pp directly --- d is always p^j for j in [0, a].
        # Recover j by counting:
        j = 0
        x = d
        while x > 1:
            x //= p
            j += 1
        sigma_d = sigma_minus_one_pp(p, j)
        convolution_value = convolution_value + mu * sigma_d**k_mpf
    assert mp.almosteq(
        formula_value,
        convolution_value,
        rel_eps=mpf("1e-30"),
        abs_eps=mpf("1e-30"),
    )


# Property 2: Factorization.from_dict is canonical ---------------------------


@given(
    pairs=st.lists(
        st.tuples(small_primes, small_exponents),
        min_size=1,
        max_size=5,
        unique_by=lambda t: t[0],
    )
)
def test_factorization_canonical(pairs):
    """Equal contents -> equal Factorization, regardless of dict order."""
    forward = Factorization.from_dict(dict(pairs))
    reverse = Factorization.from_dict(dict(reversed(pairs)))
    assert forward == reverse
    assert hash(forward) == hash(reverse)


# Property 3: Factorization.rho_k matches the pointwise rho_k -----------------


@given(
    pairs=st.lists(
        st.tuples(small_primes, small_exponents),
        min_size=1,
        max_size=4,
        unique_by=lambda t: t[0],
    )
)
def test_factorization_rho_k_matches_pointwise(pairs):
    """Computing rho^[k] from a Factorization equals pointwise rho_k(n, k)."""
    n = 1
    for p, a in pairs:
        n *= p**a
    fac = Factorization.from_dict(dict(pairs))
    k = mpf(2)
    structural = fac.rho_k(k)
    pointwise = rho_k_pointwise(n, k)
    assert mp.almosteq(
        structural,
        pointwise,
        rel_eps=mpf("1e-30"),
        abs_eps=mpf("1e-30"),
    )
