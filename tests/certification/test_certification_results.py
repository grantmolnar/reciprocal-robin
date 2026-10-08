"""Regression checks for the published bounded-parameter frontier record."""

from __future__ import annotations

import json
from decimal import Decimal
from tests._repo import REPO_ROOT

from mpmath import mpf

RESULTS = REPO_ROOT / "data" / "CERTIFICATION_RESULTS.json"


def test_advertised_frontier_record() -> None:
    data = json.loads(RESULTS.read_text(encoding="utf-8"))
    assert data["schema_version"] == "certification-results-1.0.0"
    assert data["mathematical_release"] == "reciprocal-robin-1.0.0"
    assert data["advertised_kappa"] == 100000
    run = data["frontier_run"]
    assert run["status"] == "verified"
    assert run["kappa"] == 100000
    assert run["dps"] == 60
    assert run["certificate_bytes"] == 3_624_666_835
    assert run["certificate_sha256"] == (
        "d94f3499a46c155fbcdd7e0a6f20f773ca8dd031f2af9b4217246aff9f5456c3"
    )
    assert run["tail_start_prime"] == 227977
    assert run["last_bad_prime"] == 227951
    assert run["prime_rows_checked"] == 21758
    assert run["canonical_events"] == 39454
    assert run["first_endpoint_above_5040"] == "55440"
    assert run["initial_integer_checks"] == 50399
    assert mpf(run["minimum_endpoint_log_margin_lo"]) > mpf("0.01151")
    assert run["initial_minimum_at_n"] == 10080
    assert mpf(run["initial_minimum_margin_lo"]) > mpf("0.05748")


def test_frontier_tail_summary() -> None:
    run = json.loads(RESULTS.read_text(encoding="utf-8"))["frontier_run"]
    tail = run["finite_tail_detail"]
    assert tail["A_K_floor"] == 75002
    assert tail["analytic_cutoff"] == 246590
    assert tail["tail_start_prime"] == 227977
    assert Decimal(tail["C_K_lo"]) < Decimal(tail["C_K_hi"])
    entry = run["entry_parameter"]
    assert entry["prime"] == 227977
    assert entry["exponent"] == 1
    assert Decimal(entry["epsilon_lo"]) < Decimal(entry["epsilon_hi"])
