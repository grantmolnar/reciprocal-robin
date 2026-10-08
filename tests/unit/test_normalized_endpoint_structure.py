"""Regression tests for the canonical normalized endpoint at 10800."""

from mpmath import mp, workdps

from reciprocal_robin import MpfField, using_field
from reciprocal_robin.normalized_ca import normalized_local_objective, rho_pp


def _largest_local_maximizer(p: int, eps, max_exponent: int = 12) -> int:
    values = [normalized_local_objective(p, a, 5, eps) for a in range(max_exponent + 1)]
    maximum = max(values)
    tolerance = mp.mpf("1e-55")
    return max(a for a, value in enumerate(values) if abs(value - maximum) < tolerance)


def test_10800_is_the_larger_endpoint_at_the_exact_three_adic_tie():
    with workdps(80), using_field(MpfField()):
        eps0 = mp.log(rho_pp(3, 3, 5) / rho_pp(3, 2, 5)) / mp.log(3)

        lhs = normalized_local_objective(3, 3, 5, eps0)
        rhs = normalized_local_objective(3, 2, 5, eps0)
        assert abs(lhs - rhs) < mp.mpf("1e-70")

        exponents = {p: _largest_local_maximizer(p, eps0) for p in (2, 3, 5, 7)}
        assert exponents == {2: 4, 3: 3, 5: 2, 7: 0}
        assert 2**4 * 3**3 * 5**2 == 10800


def test_local_maximizers_are_monotone_across_the_10800_transition():
    with workdps(80), using_field(MpfField()):
        eps0 = mp.log(rho_pp(3, 3, 5) / rho_pp(3, 2, 5)) / mp.log(3)
        delta = mp.mpf("1e-5")

        above = {p: _largest_local_maximizer(p, eps0 + delta) for p in (2, 3, 5, 7)}
        at = {p: _largest_local_maximizer(p, eps0) for p in (2, 3, 5, 7)}
        below = {p: _largest_local_maximizer(p, eps0 - delta) for p in (2, 3, 5, 7)}

        assert above == {2: 4, 3: 2, 5: 2, 7: 0}
        assert at == {2: 4, 3: 3, 5: 2, 7: 0}
        assert below == at
        for p in above:
            assert above[p] <= at[p] <= below[p]
