"""Tests for the retained reciprocal-CA structure theorem and safe API boundary."""

from __future__ import annotations

import pytest
from mpmath import log, mpf
from sympy import nextprime

from reciprocal_robin import F_local, MpfField, MpiField, using_field
from reciprocal_robin.exploration import (
    ThresholdTooLarge,
    audit_ca_grid,
    reciprocal_ca,
    x_threshold,
)


@pytest.mark.parametrize("k,eps", [("2", "-0.40"), ("3", "-0.25"), ("2", "0.05")])
def test_constructed_factorization_satisfies_local_bracketing(k, eps):
    kk = mpf(k)
    ee = mpf(eps)
    with using_field(MpfField()):
        fac = reciprocal_ca(kk, ee)
        for p, a in fac.primes_to_exponents.items():
            assert F_local(p, a + 1, kk) <= ee <= F_local(p, a, kk)

        # The next prime is absent, so its first increment is not profitable.
        if fac.p_max:
            q = int(nextprime(fac.p_max))
            assert F_local(q, 1, kk) <= ee


def test_local_objective_is_maximized_at_constructed_exponent():
    k = mpf(2)
    eps = mpf("-0.40")
    with using_field(MpfField()):
        fac = reciprocal_ca(k, eps)
        for p, a in fac.primes_to_exponents.items():

            def local_objective(exponent):
                from reciprocal_robin.arithmetic import rksigma_pp

                if exponent == 0:
                    return mpf(0)
                return log(rksigma_pp(p, exponent, k)) - k * eps * exponent * log(p)

            center = local_objective(a)
            assert center >= local_objective(a - 1)
            assert center >= local_objective(a + 1)


def test_exact_threshold_tie_selects_larger_local_exponent():
    with using_field(MpfField()):
        k = mpf(2)
        eps = F_local(2, 2, k)
        fac = reciprocal_ca(k, eps)
        assert fac.primes_to_exponents[2] == 2


def test_branch_dependent_ca_constructor_requires_point_field():
    with using_field(MpiField()):
        with pytest.raises(TypeError, match="point constructor"):
            reciprocal_ca(2, "-0.40")


def test_interval_grid_audit_constructs_points_then_evaluates_intervals():
    report = audit_ca_grid(2, ["-0.40"], base_dps=60, bumped_dps=120)
    assert report.all_steps_resolved
    assert len(report.steps) == 1
    assert report.steps[0].verdict in {"YES", "NO"}


def test_threshold_overflow_is_explicit_not_a_finite_sentinel():
    with using_field(MpfField()):
        with pytest.raises(ThresholdTooLarge):
            x_threshold(2, 1, mpf("-0.499999"), max_threshold=mpf(8))


def test_walk_ca_normalizes_epsilon_to_the_point_value_used():
    from mpmath import mpf
    from reciprocal_robin.exploration import walk_ca

    with using_field(MpfField()):
        step = next(walk_ca("2", ["-0.40"]))
    assert isinstance(step.eps, type(mpf("0")))
    assert step.eps == mpf("-0.40")
