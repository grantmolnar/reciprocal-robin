"""Regression tests for the retained auxiliary-function theorems."""

from __future__ import annotations

from mpmath import mpf

from reciprocal_robin import F_local, MpfField, using_field


def _strictly_decreasing(values):
    return all(left > right for left, right in zip(values, values[1:]))


def _strictly_increasing(values):
    return all(left < right for left, right in zip(values, values[1:]))


def test_F_strictly_decreases_in_x_for_real_a():
    xs = [mpf("1.01"), mpf("1.1"), mpf("1.5"), mpf(2), mpf(5), mpf(50)]
    with using_field(MpfField()):
        for k in [mpf("1.2"), mpf(2), mpf(7)]:
            for a in [mpf(1), mpf("1.5"), mpf(2), mpf(5)]:
                values = [F_local(x, a, k) for x in xs]
                assert _strictly_decreasing(values), (k, a, values)


def test_F_strictly_decreases_in_real_exponent_a():
    exponents = [mpf(1), mpf("1.25"), mpf("1.5"), mpf(2), mpf(4), mpf(10)]
    with using_field(MpfField()):
        for k in [mpf("1.2"), mpf(2), mpf(10)]:
            for x in [mpf("1.05"), mpf("1.5"), mpf(2), mpf(7)]:
                values = [F_local(x, a, k) for a in exponents]
                assert _strictly_decreasing(values), (k, x, values)


def test_F_strictly_increases_in_kappa():
    kappas = [mpf("1.01"), mpf("1.2"), mpf(2), mpf(5), mpf(20)]
    with using_field(MpfField()):
        for x in [mpf("1.1"), mpf(2), mpf(11)]:
            for a in [mpf(1), mpf("1.5"), mpf(3)]:
                values = [F_local(x, a, k) for k in kappas]
                assert _strictly_increasing(values), (x, a, values)


def test_F_limits_toward_minus_one_over_kappa():
    with using_field(MpfField()):
        for k in [mpf(2), mpf(5)]:
            target = -1 / k
            # The a=1 branch has f ~ k/x rather than 1/x, but the same F limit.
            assert abs(F_local(mpf("1e40"), 1, k) - target) < mpf("0.03")
            assert abs(F_local(mpf("1e20"), 2, k) - target) < mpf("0.03")
            assert abs(F_local(2, 100, k) - target) < mpf("1e-20")
