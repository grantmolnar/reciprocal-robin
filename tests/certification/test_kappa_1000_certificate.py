"""Regression checks for the checked-in kappa=1000 endpoint certificate."""

from __future__ import annotations

import json
from fractions import Fraction
from tests._repo import REPO_ROOT

from mpmath import exp, log, mpf


def _certificate() -> dict:
    path = REPO_ROOT / "data" / "kappa_1000_endpoint_certificate.json"
    return json.loads(path.read_text(encoding="utf-8"))


def test_kappa_1000_certificate_summary() -> None:
    data = _certificate()
    assert data["schema_version"] == "bounded-parameter-endpoints-1.0.0"
    assert data["kappa"] == "1000"

    tail = data["finite_tail"]
    assert tail["analytic_cutoff"] == 58_271
    assert tail["last_bad_prime"] == 2621
    assert tail["tail_start_prime"] == 2633
    assert tail["prime_rows_checked"] == 5902
    assert mpf(tail["terminal_block_min_margin_lo"]) > mpf("0.00029")

    endpoints = data["endpoint_verification"]
    assert endpoints["event_count"] == 794
    assert endpoints["first_endpoint_above_5040"] == "55440"
    assert endpoints["minimum_log_margin_endpoint"] == "55440"
    assert mpf(endpoints["minimum_log_margin_lo"]) > mpf("0.03909")

    initial = data["initial_integer_scan"]
    assert initial["lower"] == 5041
    assert initial["upper"] == 55439
    assert initial["count"] == 50399
    assert initial["minimum_at_n"] == 10080
    assert mpf(initial["minimum_margin_lo"]) > mpf("0.12210")


def test_kappa_1000_endpoint_prefix_starts_as_claimed() -> None:
    rows = _certificate()["endpoint_verification"]["rows"]
    endpoints = [int(row["endpoint"]) for row in rows[:9]]
    assert endpoints == [2, 6, 12, 60, 120, 360, 2520, 5040, 55440]
    assert sum(int(row["endpoint"]) > 5040 for row in rows) == 786


def test_bounded_parameter_verifier_does_not_share_package_arithmetic() -> None:
    verifier = (
        REPO_ROOT / "src" / "reciprocal_robin" / "certification" / "bounded_verify.py"
    ).read_text(encoding="utf-8")
    assert "reciprocal_robin.arithmetic" not in verifier
    assert "reciprocal_robin.normalized_ca" not in verifier
    assert "..arithmetic" not in verifier
    assert "..normalized_ca" not in verifier
    assert "bounded_generate" not in verifier


def test_kappa_1000_analytic_handoff_constants() -> None:
    assert 753 * log(753) < 5000
    assert (58_271 - 1) / 3 > 19_421
    assert exp(mpf("0.010952")) < mpf(253) / 250
    assert Fraction(71, 54) - Fraction(253, 250) == Fraction(1022, 3375)
    assert Fraction(1022, 3375) * 58_271 > 5001
