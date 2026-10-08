"""Differential tests for the incremental ratio walker.

The same CA walk is evaluated by direct recomputation and by incremental
updates. The two ratios must agree at every step.
"""

from itertools import islice

import pytest
from mpmath import mp, mpf

from reciprocal_robin import Factorization, MpfField, using_field
from reciprocal_robin.exploration import (
    CASuccessorEnumerator,
    IncrementalRatioWalker,
    SuccessorInvariantError,
)


@pytest.mark.parametrize("k", [2, 3, 5, 8, 46, 100, 1000])
def test_incremental_ratio_matches_direct(k):
    """The incremental ratio matches full recomputation at every CA step."""
    k_mpf = mpf(k)
    with using_field(MpfField()):
        # Direct recomputation.
        enum_a = CASuccessorEnumerator(k_mpf)
        direct = []
        while len(direct) < 30:
            _, fac = next(enum_a)
            if fac.log_n() <= 1:
                continue
            direct.append((fac.omega, float(fac.robin_ratio(k_mpf))))

        # Incremental update.
        walker = IncrementalRatioWalker(CASuccessorEnumerator(k_mpf), k_mpf)
        incremental = []
        for _, fac, ratio in islice(walker, 30):
            incremental.append((fac.omega, float(ratio)))

    assert len(direct) == len(incremental)
    for i, ((od, rd), (oi, ri)) in enumerate(zip(direct, incremental)):
        assert od == oi, f"step {i}: omega mismatch {od} vs {oi}"
        assert abs(rd - ri) < 1e-12 * max(abs(rd), 1), (
            f"step {i}: ratio mismatch direct={rd}, incremental={ri}, " f"diff={abs(rd - ri):.3e}"
        )


def test_incremental_walker_yields_factorization_unchanged():
    """The walker does not alter the factorizations supplied by the enumerator."""
    k = mpf(2)
    with using_field(MpfField()):
        enum = CASuccessorEnumerator(k)
        direct_steps = ((eps, fac) for eps, fac in enum if fac.log_n() > 1)
        walker = IncrementalRatioWalker(CASuccessorEnumerator(k), k)
        for (eps_a, fac_a), (eps_w, fac_w, _) in zip(islice(direct_steps, 10), islice(walker, 10)):
            assert fac_a.primes_to_exponents == fac_w.primes_to_exponents
            assert eps_a == eps_w


def test_incremental_walker_first_step():
    """The first ratio-bearing step agrees with direct recomputation."""
    k = mpf(2)
    with using_field(MpfField()):
        enum = CASuccessorEnumerator(k)
        while True:
            _, fac_direct = next(enum)
            if fac_direct.log_n() > 1:
                break
        ratio_direct = float(fac_direct.robin_ratio(k))

        walker = IncrementalRatioWalker(CASuccessorEnumerator(k), k)
        _, fac_walker, ratio_walker = next(iter(walker))

        assert fac_direct.primes_to_exponents == fac_walker.primes_to_exponents
        assert abs(ratio_direct - float(ratio_walker)) < 1e-12


def _scripted_walker(*factorizations: dict[int, int]) -> IncrementalRatioWalker:
    steps = [
        (mpf(index + 1), Factorization.from_dict(factorization))
        for index, factorization in enumerate(factorizations)
    ]
    return IncrementalRatioWalker(steps, mpf(2))


@pytest.mark.parametrize(
    ("bad_factorization", "message"),
    [
        ({2: 1, 3: 1}, "No prime changed"),
        ({2: 3, 3: 1}, "advanced by more than one exponent"),
        ({2: 2}, "Prime 3 disappeared"),
    ],
)
def test_incremental_walker_rejects_invalid_successor_steps(bad_factorization, message):
    with using_field(MpfField()):
        walker = _scripted_walker(
            {2: 1},
            {2: 1, 3: 1},
            bad_factorization,
        )
        next(walker)
        with pytest.raises(SuccessorInvariantError, match=message):
            next(walker)


def test_incremental_walker_rejects_decreasing_exponent():
    with using_field(MpfField()):
        walker = _scripted_walker(
            {2: 1},
            {2: 1, 3: 1},
            {2: 2, 3: 1},
            {2: 1, 3: 1},
        )
        next(walker)
        next(walker)
        with pytest.raises(SuccessorInvariantError, match="exponent decreased"):
            next(walker)


def test_incremental_walker_state_is_an_immutable_snapshot():
    from dataclasses import FrozenInstanceError

    with using_field(MpfField()):
        walker = _scripted_walker({2: 1}, {2: 1, 3: 1})
        _, fac, _ = next(walker)
        state = walker.state
        direct_log_n = fac.log_n()

    assert state.exponents == tuple(sorted(fac.primes_to_exponents.items()))
    assert mp.almosteq(state.log_n, direct_log_n, rel_eps=mpf("1e-30"), abs_eps=mpf("1e-30"))
    with pytest.raises(FrozenInstanceError):
        state.log_n = mpf(0)


def test_incremental_walker_rejects_backend_change_midstream():
    from reciprocal_robin import MpiField

    with using_field(MpfField()):
        walker = IncrementalRatioWalker(CASuccessorEnumerator(mpf(2)), mpf(2))
        next(walker)
        with using_field(MpiField()):
            with pytest.raises(RuntimeError, match="cannot change numerical backend"):
                next(walker)
