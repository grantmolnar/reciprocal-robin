"""Property tests for CA enumeration and the incremental ratio walker.

For random integer kappa and bounded walks, these tests check the monotone event
structure and compare the walker's running product, log n, and Robin ratio with
direct recomputation.
"""

from itertools import islice

from hypothesis import given, settings, strategies as st
from mpmath import mp, mpf

from reciprocal_robin import MpfField, using_field
from reciprocal_robin.arithmetic import rksigma_pp
from reciprocal_robin.exploration import CASuccessorEnumerator, IncrementalRatioWalker


# Integer kappa samples threshold crossings; bounded walks keep each example cheap.
k_strategy = st.integers(min_value=2, max_value=500)
steps_strategy = st.integers(min_value=1, max_value=100)


@given(k=k_strategy, n_steps=steps_strategy)
@settings(max_examples=50, deadline=None)
def test_property_eps_strictly_decreasing(k, n_steps):
    """Successive CA events have strictly decreasing eps values."""
    with using_field(MpfField()):
        enum = CASuccessorEnumerator(mpf(k))
        eps_values = [eps for eps, _ in islice(enum, n_steps)]
    for a, b in zip(eps_values, eps_values[1:]):
        assert b < a, f"k={k}: eps not decreasing: {a} -> {b}"


@given(k=k_strategy, n_steps=steps_strategy)
@settings(max_examples=50, deadline=None)
def test_property_omega_non_decreasing(k, n_steps):
    """The number of prime factors never decreases along the walk."""
    with using_field(MpfField()):
        enum = CASuccessorEnumerator(mpf(k))
        omegas = [fac.omega for _, fac in islice(enum, n_steps)]
    for a, b in zip(omegas, omegas[1:]):
        assert b >= a, f"k={k}: omega decreased: {a} -> {b}"


@given(k=k_strategy, n_steps=steps_strategy)
@settings(max_examples=50, deadline=None)
def test_property_primes_enter_in_order(k, n_steps):
    """New primes enter in increasing order."""
    with using_field(MpfField()):
        enum = CASuccessorEnumerator(mpf(k))
        seen = set()
        new_in_order = []
        for _, fac in islice(enum, n_steps):
            for p in fac.primes_to_exponents:
                if p not in seen:
                    new_in_order.append(p)
                    seen.add(p)
    assert new_in_order == sorted(
        new_in_order
    ), f"k={k}: primes did not enter in order: {new_in_order}"


@given(k=k_strategy, n_steps=steps_strategy)
@settings(max_examples=50, deadline=None)
def test_property_log_n_strictly_increasing(k, n_steps):
    """Successive CA numbers strictly increase."""
    with using_field(MpfField()):
        enum = CASuccessorEnumerator(mpf(k))
        log_ns = [fac.log_n() for _, fac in islice(enum, n_steps)]
    for a, b in zip(log_ns, log_ns[1:]):
        assert b > a, f"k={k}: log N not increasing: {a} -> {b}"


@given(k=k_strategy, n_steps=steps_strategy)
@settings(max_examples=50, deadline=None)
def test_property_per_prime_exponents_non_decreasing(k, n_steps):
    """Each prime exponent is nondecreasing along the walk."""
    with using_field(MpfField()):
        enum = CASuccessorEnumerator(mpf(k))
        factorizations = [fac.primes_to_exponents for _, fac in islice(enum, n_steps)]
    all_primes = set()
    for f in factorizations:
        all_primes.update(f.keys())
    for p in all_primes:
        # Absence means exponent zero.
        seq = [f.get(p, 0) for f in factorizations]
        for a, b in zip(seq, seq[1:]):
            assert b >= a, f"k={k}: prime {p}'s exponent decreased: {a} -> {b}"


@given(k=k_strategy, n_steps=steps_strategy)
@settings(max_examples=50, deadline=None)
def test_property_walker_prod_matches_direct_formula(k, n_steps):
    """The stored product equals the direct product at every step."""
    k_mpf = mpf(k)
    with using_field(MpfField()):
        walker = IncrementalRatioWalker(CASuccessorEnumerator(k_mpf), k_mpf)
        for step_idx, (_, fac, _) in enumerate(islice(walker, n_steps)):
            # Recompute the product from the current factorization.
            direct_prod = mpf(1)
            for p, a in fac.primes_to_exponents.items():
                direct_prod = direct_prod * mpf(p) ** a * rksigma_pp(p, a, k_mpf) / k_mpf
            incremental_prod = walker.state.product
            scale = max(abs(direct_prod), mpf("1e-50"))
            rel_diff = abs(direct_prod - incremental_prod) / scale
            assert rel_diff < mpf("1e-30"), (
                f"k={k}, step {step_idx}: prod mismatch, "
                f"direct={mp.nstr(direct_prod, 12)}, "
                f"incremental={mp.nstr(incremental_prod, 12)}, "
                f"rel diff={mp.nstr(rel_diff, 6)}"
            )


@given(k=k_strategy, n_steps=steps_strategy)
@settings(max_examples=50, deadline=None)
def test_property_walker_log_n_matches_fac_log_n(k, n_steps):
    """The stored log n equals direct recomputation at every step."""
    k_mpf = mpf(k)
    with using_field(MpfField()):
        walker = IncrementalRatioWalker(CASuccessorEnumerator(k_mpf), k_mpf)
        for step_idx, (_, fac, _) in enumerate(islice(walker, n_steps)):
            direct_log_n = fac.log_n()
            incremental_log_n = walker.state.log_n
            assert mp.almosteq(
                direct_log_n,
                incremental_log_n,
                rel_eps=mpf("1e-30"),
                abs_eps=mpf("1e-30"),
            ), (
                f"k={k}, step {step_idx}: log_n mismatch, "
                f"direct={mp.nstr(direct_log_n, 12)}, "
                f"incremental={mp.nstr(incremental_log_n, 12)}"
            )


@given(k=k_strategy, n_steps=steps_strategy)
@settings(max_examples=50, deadline=None)
def test_property_walker_ratio_matches_fac_robin_ratio(k, n_steps):
    """The incremental Robin ratio agrees with direct recomputation."""
    k_mpf = mpf(k)
    with using_field(MpfField()):
        walker = IncrementalRatioWalker(CASuccessorEnumerator(k_mpf), k_mpf)
        for step_idx, (_, fac, ratio_w) in enumerate(islice(walker, n_steps)):
            ratio_d = fac.robin_ratio(k_mpf)
            scale = max(abs(ratio_d), mpf("1e-50"))
            rel_diff = abs(ratio_d - ratio_w) / scale
            assert rel_diff < mpf("1e-20"), (
                f"k={k}, step {step_idx}: ratio mismatch, "
                f"direct={mp.nstr(ratio_d, 12)}, "
                f"incremental={mp.nstr(ratio_w, 12)}, "
                f"rel diff={mp.nstr(rel_diff, 6)}"
            )
