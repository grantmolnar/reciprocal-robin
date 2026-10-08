from pathlib import Path

from tests._repo import REPO_ROOT
import json

import pytest

from mpmath import mpf

from reciprocal_robin import MpfField, using_field
from reciprocal_robin.normalized_ca import rho_pp, square_entry_gap_k5
from reciprocal_robin._serialization import write_json
from reciprocal_robin.certification.kappa_five_generate import main as generate_kappa_five
from reciprocal_robin.certification.kappa_five_verify import verify_certificate


def test_square_entry_gap_positive_for_sample_primes():
    with using_field(MpfField()):
        for p in (2, 3, 5, 7, 11, 101, 1009):
            assert square_entry_gap_k5(p) > 0
            assert rho_pp(p, 2, 5) > rho_pp(p, 1, 5) ** 2


def test_square_entry_changes_direction_above_five_at_large_prime():
    with using_field(MpfField()):
        p = 1009
        assert rho_pp(p, 2, mpf(6)) < rho_pp(p, 1, mpf(6)) ** 2


def test_checked_in_certificate_passes():
    path = REPO_ROOT / "data" / "kappa_five_interval_certificate.json"
    data = json.loads(path.read_text())
    assert data["mertens_margin"]["all_positive"] is True
    assert data["mertens_margin"]["minimum_at_prime"] == 11
    threshold = data["primorial_threshold"]
    assert threshold["certified_positive"] is True
    assert mpf(threshold["margin_lo"]) > mpf("0.65")
    extension = data["mertens_extension"]
    assert extension["correction_constant"] == "0.10952"
    assert extension["prime_count"] == 165866
    assert extension["all_verified"] is True
    assert extension["minimum_log_margin_at_prime"] == 24137
    assert extension["maximum_required_correction_at_prime"] == 24137
    required = mpf(extension["maximum_required_correction_upper"])
    assert required > mpf("0.10836")
    assert required < mpf("0.10952")
    assert data["normalized_ca_seed"]["largest_endpoint"] == 10800
    assert data["normalized_ca_seed"]["all_checks_pass"] is True
    checks = data["normalized_ca_seed"]["checks"]
    assert checks["p3_a2_a3_tie_exact_by_definition"] is True
    assert checks["p3_a2_a3_interval_consistency"] is True
    assert data["low_support_margin"]["certified_positive"] is True


def test_kappa_five_verifier_is_independent_of_generator_arithmetic() -> None:
    verifier = (
        REPO_ROOT / "src" / "reciprocal_robin" / "certification" / "kappa_five_verify.py"
    ).read_text(encoding="utf-8")
    assert "kappa_five_generate" not in verifier
    assert "..arithmetic" not in verifier
    assert "..normalized_ca" not in verifier


def test_kappa_five_verifier_rejects_malformed_header_before_reconstruction(
    tmp_path: Path,
) -> None:
    certificate = {
        "schema_version": "kappa-five-interval-1.1.0",
        "status": "directed-interval-finite-certificate",
        "precision_dps": 60,
        "theorem_scope": (
            "shared finite inputs for the unconditional interval 1 < kappa <= 5 "
            "and the uniform explicit upper bound"
        ),
    }
    mutations = [
        lambda data: data.__setitem__("schema_version", "wrong-schema"),
        lambda data: data.__setitem__("status", "not-a-certificate"),
        lambda data: data.__setitem__("theorem_scope", "wrong theorem"),
        lambda data: data.__setitem__("precision_dps", "60"),
        lambda data: data.__setitem__("precision_dps", 60.0),
        lambda data: data.__setitem__("precision_dps", True),
        lambda data: data.__setitem__("precision_dps", 1),
    ]
    for index, mutate in enumerate(mutations):
        tampered = certificate.copy()
        mutate(tampered)
        path = tmp_path / f"bad-header-{index}.json"
        write_json(path, tampered)
        with pytest.raises(AssertionError):
            verify_certificate(path)

    non_object = tmp_path / "non-object.json"
    write_json(non_object, [])
    with pytest.raises(AssertionError, match="JSON object"):
        verify_certificate(non_object)


def test_kappa_five_generator_rejects_invalid_cli_precision() -> None:
    for dps in ("1", "0", "-10", "not-an-integer"):
        with pytest.raises(SystemExit) as excinfo:
            generate_kappa_five(["--dps", dps])
        assert excinfo.value.code == 2
