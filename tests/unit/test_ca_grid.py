"""Tests for sampled reciprocal-CA grid reports."""

import json

from mpmath import mpf

from reciprocal_robin.exploration import CAWalkReport, audit_ca_grid

EPS_GRID = ["-0.30", "-0.35", "-0.40", "-0.43", "-0.45", "-0.47"]


def test_ca_grid_report_resolves_sampled_points():
    report = audit_ca_grid(mpf(2), [mpf(e) for e in EPS_GRID])
    assert isinstance(report, CAWalkReport)
    assert report.all_steps_resolved is True
    assert report.ambiguity_witness is None
    assert len(report.steps) == len(EPS_GRID)
    assert all(step.verdict == "YES" for step in report.steps)


def test_ca_grid_report_records_corrected_structure():
    report = audit_ca_grid(mpf(2), [mpf("-0.40")])
    step = report.steps[0]
    assert step.omega == 11
    assert step.p_max == 31
    assert float(step.ratio_hi) < 0.463


def test_ca_grid_report_json_schema_has_no_universal_bounds():
    report = audit_ca_grid(mpf(2), [mpf("-0.40")])
    parsed = json.loads(report.to_json())
    assert parsed["schema_version"] == "ca-walk-report-1.1.0"
    assert parsed["all_steps_resolved"] is True
    assert "n_lower_bound" not in parsed
    assert "log_n" not in parsed["steps"][0]
    assert float(parsed["steps"][0]["log_n_lo"]) <= float(parsed["steps"][0]["log_n_hi"])


def test_report_intervals_are_well_formed():
    report = audit_ca_grid(mpf(2), [mpf(e) for e in EPS_GRID])
    for step in report.steps:
        assert float(step.ratio_lo) <= float(step.ratio_hi)


def test_ca_grid_report_steps_are_immutable():
    report = audit_ca_grid("2", ["-0.40"])
    assert isinstance(report.steps, tuple)


def test_text_inputs_are_parsed_under_requested_precision():
    from mpmath import mp

    prior = mp.dps
    try:
        mp.dps = 15
        report = audit_ca_grid(
            "2.0000000000000000000000000000000000000001",
            ["-0.4000000000000000000000000000000000000001"],
            base_dps=60,
            bumped_dps=120,
        )
    finally:
        mp.dps = prior
    assert report.kappa.startswith("2.0000000000000000000000000000000000000001")


def test_grid_helpers_reject_empty_or_vacuous_samples():
    import pytest
    from reciprocal_robin import MpfField, using_field
    from reciprocal_robin.exploration import max_robin_ratio_on_ca_grid

    with using_field(MpfField()):
        with pytest.raises(ValueError, match="at least one"):
            max_robin_ratio_on_ca_grid(2, [])
        with pytest.raises(ValueError, match="no nontrivial"):
            max_robin_ratio_on_ca_grid(2, ["100"])
    with pytest.raises(ValueError, match="at least one"):
        audit_ca_grid(2, [])
    with pytest.raises(ValueError, match="no nontrivial"):
        audit_ca_grid(2, ["100"])
