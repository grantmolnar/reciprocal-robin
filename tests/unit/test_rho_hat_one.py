"""Tests for the boundary limit rho^[1](n) = lim_{k->1+} rho^[k](n).

These tests pin ``MpfField`` because they benchmark computed values against
known reference numbers via point comparisons. A separate test below
verifies that ``rho_hat_one`` runs correctly under the default ``MpiField``
and produces intervals that bracket the point value.
"""

from mpmath import mpf, euler, log, exp

from reciprocal_robin import rho_hat_one, rho_k, MpfField, MpiField, using_field


def test_approaches_limit_from_above():
    """rho^[k](n) should approach rho^[1](n) as k -> 1+, monotonically from above."""
    with using_field(MpfField()):
        for n in [3, 6, 12, 24, 36, 5040]:
            limit = rho_hat_one(n)
            vals = [rho_k(n, mpf(f"1.{'0'*j}1")) for j in [2, 3, 4, 5, 6]]
            for i in range(len(vals) - 1):
                assert (
                    vals[i] > vals[i + 1] > limit
                ), f"n={n}: rho at k_i={vals[i]}, k_{i+1}={vals[i+1]}, limit={limit}"
            assert abs(vals[-1] - limit) / abs(limit) < 1e-3, f"n={n}: {vals[-1]} vs {limit}"


def test_strictly_less_than_sigma_minus_one():
    """Fundamental lemma: rho^[k](n) < sigma_{-1}(n) for all k >= 1, all n > 1."""
    from sympy import divisors

    with using_field(MpfField()):
        for n in [3, 6, 12, 24, 36, 100, 5040]:
            sigma_m1 = sum(mpf(1) / mpf(d) for d in divisors(n))
            assert rho_hat_one(n) < sigma_m1, f"n={n}: rho^[1] >= sigma_{{-1}}"


def test_boundary_violator_characterization():
    """rho^[1](n) >= e^gamma log log n iff n in {3,4,5,6,8,12,16,24}."""
    expected = {3, 4, 5, 6, 8, 12, 16, 24}
    candidates = [
        3,
        4,
        5,
        6,
        8,
        9,
        10,
        12,
        16,
        18,
        20,
        24,
        30,
        36,
        48,
        60,
        72,
        84,
        120,
        180,
        240,
        360,
        720,
        840,
        2520,
        5040,
    ]
    with using_field(MpfField()):
        boundary = {n for n in candidates if rho_hat_one(n) >= exp(euler) * log(log(n))}
    assert boundary == expected, f"got {sorted(boundary)}, expected {sorted(expected)}"


def test_multiplicativity():
    """rho^[1] is multiplicative on coprime arguments."""
    from mpmath import almosteq

    with using_field(MpfField()):
        lhs = rho_hat_one(8 * 9)
        rhs = rho_hat_one(8) * rho_hat_one(9)
        assert almosteq(lhs, rhs, rel_eps=mpf("1e-30")), f"{lhs} vs {rhs}"

        lhs = rho_hat_one(35)
        rhs = rho_hat_one(5) * rho_hat_one(7)
        assert almosteq(lhs, rhs, rel_eps=mpf("1e-30"))


def test_rho_hat_one_under_mpi_field_produces_interval_bracketing_point_value():
    """Under MpiField, rho_hat_one returns an interval bracketing the point value."""
    import mpmath

    with using_field(MpfField()):
        point_values = {n: float(rho_hat_one(n)) for n in [3, 12, 24, 36, 5040]}
    with using_field(MpiField()):
        for n, expected in point_values.items():
            interval = rho_hat_one(n)
            lo = float(mpmath.mpf(interval.a))
            hi = float(mpmath.mpf(interval.b))
            assert lo <= expected <= hi, f"n={n}: interval [{lo}, {hi}] does not bracket {expected}"
