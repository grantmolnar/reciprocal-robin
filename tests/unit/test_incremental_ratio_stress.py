"""Stress tests for the incremental ratio walker.

Long walks test accumulated point error, interval widening, and possible
collisions among local successor parameters.
"""

from itertools import islice

import pytest
from mpmath import mpf

from reciprocal_robin import MpfField, MpiField, using_field
from reciprocal_robin.arithmetic import F_local
from reciprocal_robin.exploration import (
    CASuccessorEnumerator,
    IncrementalRatioWalker,
    SuccessorInvariantError,
)


pytestmark = pytest.mark.stress


def _ratio_bearing_steps(enum):
    """Yield successor steps only after the structural n=2 event."""
    for eps, fac in enum:
        if fac.log_n() > 1:
            yield eps, fac


# --- Drift over long walks --------------------------------------------------


@pytest.mark.parametrize(
    "k,n_steps",
    [
        (2, 500),
        (3, 500),
        (100, 300),
    ],
)
def test_no_drift_over_long_walk(k, n_steps):
    """Long point-arithmetic walks remain close to direct recomputation."""
    k_mpf = mpf(k)
    with using_field(MpfField()):
        enum_direct = CASuccessorEnumerator(k_mpf)
        walker = IncrementalRatioWalker(CASuccessorEnumerator(k_mpf), k_mpf)

        max_rel_diff = 0.0
        steps_checked = 0
        for (_, fac_d), (_, fac_w, ratio_w) in zip(
            islice(_ratio_bearing_steps(enum_direct), n_steps), islice(walker, n_steps)
        ):
            direct = float(fac_d.robin_ratio(k_mpf))
            incremental = float(ratio_w)
            # Guard the relative error near zero.
            denom = max(abs(direct), 1e-30)
            rel = abs(direct - incremental) / denom
            max_rel_diff = max(max_rel_diff, rel)
            steps_checked += 1

        # The threshold is loose relative to 50-digit arithmetic.
        assert steps_checked == n_steps, f"Only checked {steps_checked}/{n_steps} steps"
        assert max_rel_diff < 1e-10, (
            f"Drift over {n_steps} steps at k={k}: " f"max relative diff = {max_rel_diff:.3e}"
        )


# --- Interval width under MpiField ------------------------------------------


@pytest.mark.parametrize(
    "k,n_steps",
    [
        (2, 100),
        (8, 100),
    ],
)
def test_interval_width_stays_narrow(k, n_steps):
    """Incremental intervals remain comparable in width to direct intervals."""
    k_mpi = mpf(k)
    with using_field(MpiField()):
        enum_direct = CASuccessorEnumerator(k_mpi)
        walker = IncrementalRatioWalker(CASuccessorEnumerator(k_mpi), k_mpi)

        direct_widths = []
        inc_widths = []
        for (_, fac_d), (_, fac_w, ratio_w) in zip(
            islice(_ratio_bearing_steps(enum_direct), n_steps), islice(walker, n_steps)
        ):
            direct_ratio = fac_d.robin_ratio(k_mpi)
            dw = float(direct_ratio.b) - float(direct_ratio.a)
            iw = float(ratio_w.b) - float(ratio_w.a)
            direct_widths.append(dw)
            inc_widths.append(iw)

        # Four decimal digits of additional widening is the allowed ceiling.
        ratios = [iw / dw if dw > 0 else 1.0 for dw, iw in zip(direct_widths, inc_widths)]
        max_ratio = max(ratios)
        assert max_ratio < 1e4, (
            f"Incremental interval {max_ratio:.1f}x wider than direct at k={k} "
            f"after {n_steps} steps. Direct max width = "
            f"{max(direct_widths):.3e}, incremental max width = "
            f"{max(inc_widths):.3e}"
        )


# --- F^[k] collision search -------------------------------------------------


def test_F_local_no_collisions_in_reasonable_range():
    """Search the tested range for collisions among local successor parameters.

    The check is empirical rather than a proof of non-collision. Point arithmetic
    avoids the false positives introduced by float truncation when nearby values
    approach the same limit.
    """
    import mpmath
    from sympy import primerange

    primes = list(primerange(2, 101))
    EPS_TOL = mpmath.mpf("1e-40")  # Numerical equality threshold at 50 dps.

    for k in [mpf(2), mpf(3), mpf(8), mpf(100)]:
        with using_field(MpfField()):
            values = []
            for p in primes:
                for ell in range(1, 21):
                    eps = F_local(p, ell, k)
                    values.append((eps, p, ell))
            # A collision must occur between adjacent values after sorting.
            values.sort(key=lambda t: t[0])
            for (e1, p1, l1), (e2, p2, l2) in zip(values, values[1:]):
                if abs(e1 - e2) < EPS_TOL:
                    pytest.fail(
                        f"F_local collision at k={k}: "
                        f"F^[k]({p1},{l1}) ~ F^[k]({p2},{l2}) "
                        f"(diff = {abs(e1-e2)})"
                    )


def test_walker_detects_multi_prime_change():
    """A successor stream that changes two primes violates the walker contract."""
    from reciprocal_robin import Factorization

    k = mpf(2)
    scripted_steps = iter(
        [
            (mpf("0.5"), Factorization.from_dict({2: 1})),
            (mpf("0.4"), Factorization.from_dict({2: 1, 3: 1})),
            (mpf("0.3"), Factorization.from_dict({2: 2, 3: 1, 5: 1})),
        ]
    )
    with using_field(MpfField()):
        walker = IncrementalRatioWalker(scripted_steps, k)
        next(walker)
        with pytest.raises(SuccessorInvariantError, match="More than one prime"):
            next(walker)
