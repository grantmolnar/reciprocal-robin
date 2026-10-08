"""Tests for numerical fields, precision state, and the Trilean enum."""

import os
import subprocess
import sys

import pytest
from mpmath import iv, mp, mpf

from tests._repo import REPO_ROOT

from reciprocal_robin import (
    MpfField,
    RealField,
    Trilean,
    current_field,
    robin_ratio,
    using_field,
    working_precision,
)


def test_mpf_field_satisfies_protocol():
    assert isinstance(MpfField(), RealField)


def test_default_field_is_mpi():
    """The default interval field permits an ambiguous comparison result."""
    assert current_field().name == "mpi"


def test_using_field_swaps_and_restores():
    outer = current_field()
    with using_field(MpfField()) as inner:
        assert current_field() is inner
        assert inner is not outer or current_field().name == "mpf"
    assert current_field() is outer


def test_using_field_restores_on_exception():
    outer = current_field()
    try:
        with using_field(MpfField()):
            raise RuntimeError("boom")
    except RuntimeError:
        pass
    assert current_field() is outer


def test_trilean_has_three_values():
    assert {Trilean.YES, Trilean.NO, Trilean.AMBIGUOUS} == set(Trilean)


def test_mpf_field_trilean_is_binary():
    F = MpfField()
    assert F.is_strictly_less_than_one(F.number("0.5")) == Trilean.YES
    assert F.is_strictly_less_than_one(F.number("1.5")) == Trilean.NO
    # Point arithmetic has no ambiguous comparison state.
    assert F.is_strictly_less_than_one(F.number(1)) == Trilean.NO


def test_arithmetic_through_explicit_field_matches_default():
    """Point evaluation is deterministic, and interval evaluation contains it."""
    with using_field(MpfField()):
        pinned_mpf = float(robin_ratio(36, mpf(2)))
    with using_field(MpfField()):
        pinned_mpf_again = float(robin_ratio(36, mpf(2)))
    assert pinned_mpf == pinned_mpf_again

    # The default interval field contains the point value.
    default_interval = robin_ratio(36, mpf(2))
    assert float(default_interval.a) <= pinned_mpf <= float(default_interval.b)


# ----- MpiField (interval arithmetic) ---------------------------------------


def test_mpi_field_satisfies_protocol():
    from reciprocal_robin import MpiField

    assert isinstance(MpiField(), RealField)


def test_mpi_field_trilean_three_regimes():
    """An interval clearly below 1, clearly above 1, and straddling 1."""
    from reciprocal_robin import MpiField

    F = MpiField()
    assert F.is_strictly_less_than_one(F.number("0.5")) == Trilean.YES
    assert F.is_strictly_less_than_one(F.number("2")) == Trilean.NO
    # This interval genuinely straddles 1.
    import mpmath

    straddling = mpmath.mpi("0.9", "1.1")
    assert F.is_strictly_less_than_one(straddling) == Trilean.AMBIGUOUS


def test_robin_ratio_through_mpi_field_brackets_mpf_value():
    """The interval evaluation of robin_ratio(36, 2) contains the point value."""
    from reciprocal_robin import MpiField

    with using_field(MpfField()):
        point_value = float(robin_ratio(36, mpf(2)))
    with using_field(MpiField()):
        interval = robin_ratio(36, mpf(2))
    # Interval containment is the relevant comparison.
    assert float(interval.a) <= point_value <= float(interval.b)
    # The default precision also resolves this case cleanly.
    F = MpiField()
    assert F.is_strictly_less_than_one(interval) == Trilean.YES


def test_working_precision_actually_tightens_intervals():
    """Increasing scoped precision narrows interval evaluations.

    mpmath keeps point and interval precision separately; the package context
    updates both.
    """
    from reciprocal_robin import MpiField, working_precision

    widths = []
    for dps in (30, 100):
        with working_precision(dps), using_field(MpiField()):
            r = robin_ratio(36, mpf(2))
            widths.append(float(r.b - r.a))
    assert widths[1] < widths[0], f"interval width should shrink with precision, got {widths}"


def test_import_does_not_mutate_mpmath_precision() -> None:
    code = """
from mpmath import iv, mp
mp.dps = 37
iv.dps = 41
before = (mp.dps, iv.dps)
import reciprocal_robin
assert (mp.dps, iv.dps) == before
"""
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO_ROOT / "src")
    subprocess.run([sys.executable, "-c", code], env=env, check=True)


def test_working_precision_restores_both_contexts() -> None:
    before = (mp.dps, iv.dps)
    with working_precision(73):
        assert (mp.dps, iv.dps) == (73, 73)
    assert (mp.dps, iv.dps) == before


@pytest.mark.parametrize("bad_dps", [0, 1, -5, 2.5, True])
def test_working_precision_rejects_invalid_values(bad_dps) -> None:
    with pytest.raises((TypeError, ValueError)):
        with working_precision(bad_dps):
            pass
