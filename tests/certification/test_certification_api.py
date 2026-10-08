"""Public architecture tests for bounded-parameter certification tools."""

from __future__ import annotations

import argparse
import ast
import sys
from pathlib import Path

from tests._repo import REPO_ROOT

from reciprocal_robin._serialization import write_json
from reciprocal_robin.certification.bounded_generate import compute_certificate
from reciprocal_robin.certification.bounded_verify import verify_certificate
from reciprocal_robin.certification.sweep import geometric_values, parse_values
from reciprocal_robin.certification.runtime import unbounded_decimal_integers


def _package_imports(path: Path) -> set[str]:
    """Resolve package-local imports without importing the module under test."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imports: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(
                alias.name for alias in node.names if alias.name.startswith("reciprocal_robin")
            )
            continue
        if not isinstance(node, ast.ImportFrom):
            continue
        if node.level == 0:
            if node.module and node.module.startswith("reciprocal_robin"):
                imports.add(node.module)
            continue
        if node.level == 1:
            base = "reciprocal_robin.certification"
        elif node.level == 2:
            base = "reciprocal_robin"
        else:
            continue
        imports.add(f"{base}.{node.module}" if node.module else base)
    return imports


def test_independent_verifiers_have_only_reviewed_package_dependencies() -> None:
    certification = REPO_ROOT / "src" / "reciprocal_robin" / "certification"
    bounded_imports = _package_imports(certification / "bounded_verify.py")
    kappa_five_imports = _package_imports(certification / "kappa_five_verify.py")

    assert bounded_imports == {
        "reciprocal_robin._serialization",
        "reciprocal_robin.certification.runtime",
    }
    assert kappa_five_imports == {"reciprocal_robin._serialization"}


def test_generator_and_verifier_round_trip_small_certificate(tmp_path: Path) -> None:
    certificate = tmp_path / "kappa_10.json"
    write_json(certificate, compute_certificate("10", dps=40))
    report = verify_certificate(certificate)
    assert report["all_checks_pass"] is True
    assert report["kappa"] == "10"


def test_sweep_value_helpers() -> None:
    assert parse_values("2000, 4000,8000") == [2000, 4000, 8000]
    assert list(geometric_values(2000, 2.0, 8000)) == [2000, 4000, 8000]


def test_installed_bounded_generator_requires_kappa_argument() -> None:
    import pytest

    from reciprocal_robin.certification.bounded_generate import main

    with pytest.raises(SystemExit) as excinfo:
        main([])
    assert excinfo.value.code == 2


def test_installed_bounded_verifier_requires_certificate_argument() -> None:
    import pytest

    from reciprocal_robin.certification.bounded_verify import main

    with pytest.raises(SystemExit) as excinfo:
        main([])
    assert excinfo.value.code == 2


def test_installed_kappa_five_verifier_requires_certificate_argument() -> None:
    import pytest

    from reciprocal_robin.certification.kappa_five_verify import main

    with pytest.raises(SystemExit) as excinfo:
        main([])
    assert excinfo.value.code == 2


def test_certification_runtime_allows_large_decimal_endpoints_and_restores_limit() -> None:
    getter = getattr(sys, "get_int_max_str_digits", None)
    setter = getattr(sys, "set_int_max_str_digits", None)
    if getter is None or setter is None:
        return

    previous = getter()
    try:
        setter(4300)
        endpoint = 10**5000
        with unbounded_decimal_integers():
            text = str(endpoint)
            assert len(text) == 5001
            assert int(text) == endpoint
        assert getter() == 4300
    finally:
        setter(previous)


def test_certification_runtime_restores_limit_after_nested_exception() -> None:
    getter = getattr(sys, "get_int_max_str_digits", None)
    setter = getattr(sys, "set_int_max_str_digits", None)
    if getter is None or setter is None:
        return

    import pytest

    previous = getter()
    try:
        setter(4300)
        with pytest.raises(RuntimeError, match="sentinel"):
            with unbounded_decimal_integers():
                assert getter() == 0
                with unbounded_decimal_integers():
                    assert getter() == 0
                assert getter() == 0
                raise RuntimeError("sentinel")
        assert getter() == 4300
    finally:
        setter(previous)


def test_bounded_verifier_rejects_contradictory_discrete_metadata(tmp_path: Path) -> None:
    import copy
    import pytest

    certificate = compute_certificate("10", dps=40)

    mutations = [
        ("status", lambda data: data.__setitem__("status", "not-a-certificate")),
        (
            "endpoint positivity",
            lambda data: data["endpoint_verification"].__setitem__("all_positive", False),
        ),
        (
            "initial positivity",
            lambda data: data["initial_integer_scan"].__setitem__("all_positive", False),
        ),
        (
            "endpoint event metadata",
            lambda data: data["endpoint_verification"]["rows"][0].__setitem__("event_prime", 99991),
        ),
        (
            "hull order",
            lambda data: data["local_profiles"][0].__setitem__(
                "hull_vertices", [0, 1, 1, data["local_profiles"][0]["exponent_cap"]]
            ),
        ),
    ]

    for label, mutate in mutations:
        tampered = copy.deepcopy(certificate)
        mutate(tampered)
        path = tmp_path / f"tampered-{label.replace(' ', '-')}.json"
        write_json(path, tampered)
        with pytest.raises(AssertionError):
            verify_certificate(path)


def test_unified_certification_dispatcher_does_not_eagerly_import_proof_modules() -> None:
    import os
    import subprocess

    repo_root = REPO_ROOT
    code = """
import sys
import reciprocal_robin.certification.__main__
for name in (
    'reciprocal_robin.certification.bounded_generate',
    'reciprocal_robin.certification.bounded_verify',
    'reciprocal_robin.certification.kappa_five_generate',
    'reciprocal_robin.certification.kappa_five_verify',
):
    assert name not in sys.modules, name
"""
    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo_root / "src")
    subprocess.run([sys.executable, "-c", code], env=env, check=True)


def test_sweep_rejects_invalid_control_parameters_without_running_stages() -> None:
    import pytest

    from reciprocal_robin.certification.sweep import main

    invalid = [
        ["--dps", "1"],
        ["--timeout-minutes", "-1"],
        ["--timeout-minutes", "nan"],
        ["--factor", "nan"],
        ["--factor", "1"],
    ]
    for argv in invalid:
        with pytest.raises(SystemExit) as excinfo:
            main(argv)
        assert excinfo.value.code == 2


def test_parse_values_reports_nonintegers_cleanly() -> None:
    import pytest

    with pytest.raises(argparse.ArgumentTypeError, match="must be an integer"):
        parse_values("2000,nope,8000")


def test_unified_certification_dispatcher_has_documented_command_surface() -> None:
    from reciprocal_robin.certification import __main__ as dispatcher

    assert tuple(dispatcher._COMMANDS) == (
        "generate",
        "verify",
        "sweep",
        "generate-kappa-five",
        "verify-kappa-five",
    )
