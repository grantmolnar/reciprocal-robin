"""Structural tests for the certificate and release inventory."""

from __future__ import annotations

import json
import os
from tests._repo import REPO_ROOT
import subprocess
import sys

MANIFEST = REPO_ROOT / "PROOF_MANIFEST.json"


def test_manifest_lists_existing_artifacts() -> None:
    """Check manifest structure without comparing release hashes."""
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert data["mathematical_release"] == "reciprocal-robin-1.0.0"
    assert data["package_version"] == "1.0.0"
    for item in data["artifacts"]:
        path = REPO_ROOT / item["path"]
        assert path.is_file()
        assert isinstance(item["bytes"], int) and item["bytes"] >= 0
        assert isinstance(item["sha256"], str) and len(item["sha256"]) == 64


def test_manifest_has_expected_load_bearing_data_artifacts():
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    load_bearing = [item for item in data["artifacts"] if item["load_bearing"]]
    assert [item["path"] for item in load_bearing] == [
        "data/kappa_five_interval_certificate.json",
        "data/CERTIFICATION_RESULTS.json",
    ]


def test_manifest_includes_certification_subpackage() -> None:
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    paths = {item["path"] for item in data["artifacts"]}
    assert "src/reciprocal_robin/certification/bounded_generate.py" in paths
    assert "src/reciprocal_robin/certification/bounded_verify.py" in paths
    assert "src/reciprocal_robin/certification/sweep.py" in paths
    assert {"README.md", "ARCHITECTURE.md", "CERTIFICATION.md"} <= paths


def test_package_version_is_consistent_across_release_metadata():
    import tomllib

    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        pyproject_version = tomllib.load(handle)["tool"]["poetry"]["version"]
    citation_lines = (REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8").splitlines()
    citation_version = next(
        line.split(":", 1)[1].strip() for line in citation_lines if line.startswith("version:")
    )
    assert data["package_version"] == pyproject_version == citation_version


def test_portable_verification_pins_match_poetry_lock() -> None:
    import re
    import tomllib

    with (REPO_ROOT / "poetry.lock").open("rb") as handle:
        lock = tomllib.load(handle)
    locked = {
        package["name"].lower().replace("_", "-"): package["version"] for package in lock["package"]
    }

    pins: dict[str, str] = {}
    for raw in (
        (REPO_ROOT / "requirements-verification.txt").read_text(encoding="utf-8").splitlines()
    ):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        match = re.fullmatch(r"([A-Za-z0-9_.-]+)==([^\s]+)", line)
        assert match is not None, f"verification requirement is not an exact pin: {line}"
        pins[match.group(1).lower().replace("_", "-")] = match.group(2)

    assert pins
    assert {name: locked.get(name) for name in pins} == pins


def test_manifest_matches_current_builder_inventory() -> None:
    """Check inventory membership independently of byte-level release freshness."""
    scripts = REPO_ROOT / "scripts"
    sys.path.insert(0, str(scripts))
    try:
        from build_proof_manifest import specifications
    finally:
        sys.path.remove(str(scripts))

    checked_in = json.loads(MANIFEST.read_text(encoding="utf-8"))
    checked_paths = {item["path"] for item in checked_in["artifacts"]}
    expected_paths = {item["path"] for item in specifications()}
    assert checked_paths == expected_paths


def test_release_verifier_rejects_stale_manifest(tmp_path, monkeypatch) -> None:
    scripts = REPO_ROOT / "scripts"
    sys.path.insert(0, str(scripts))
    try:
        import verify_release
    finally:
        sys.path.remove(str(scripts))

    stale = json.loads(MANIFEST.read_text(encoding="utf-8"))
    stale["artifacts"] = stale["artifacts"][:-1]
    path = tmp_path / "stale-manifest.json"
    path.write_text(json.dumps(stale), encoding="utf-8")
    monkeypatch.setattr(verify_release, "MANIFEST_PATH", path)
    report, ok = verify_release.verify_manifest()
    assert ok is False
    assert report["inventory_matches_builder"] is False


def test_release_pin_parser_rejects_malformed_and_duplicate_entries(tmp_path, monkeypatch) -> None:
    scripts = REPO_ROOT / "scripts"
    sys.path.insert(0, str(scripts))
    try:
        import verify_release
    finally:
        sys.path.remove(str(scripts))

    malformed = tmp_path / "malformed.txt"
    malformed.write_text("mpmath>=1.3.0\n", encoding="utf-8")
    monkeypatch.setattr(verify_release, "REQUIREMENTS_PATH", malformed)
    import pytest

    with pytest.raises(ValueError, match="invalid exact pin"):
        verify_release.parse_pins()

    duplicate = tmp_path / "duplicate.txt"
    duplicate.write_text("mpmath==1.3.0\nMPMATH==1.3.0\n", encoding="utf-8")
    monkeypatch.setattr(verify_release, "REQUIREMENTS_PATH", duplicate)
    with pytest.raises(ValueError, match="duplicate verification pin"):
        verify_release.parse_pins()


def test_proof_scripts_bootstrap_from_a_source_checkout() -> None:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    scripts = [
        "certify_kappa_five_interval.py",
        "verify_kappa_five_certificate.py",
        "certify_bounded_parameter_interval.py",
        "verify_bounded_parameter_certificate.py",
        "sweep_bounded_parameter.py",
        "verify_release.py",
    ]
    for script in scripts:
        proc = subprocess.run(
            [sys.executable, str(REPO_ROOT / "scripts" / script), "--help"],
            cwd=REPO_ROOT,
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        assert proc.returncode == 0, f"{script}: {proc.stderr}"


def test_manuscript_runners_bootstrap_from_a_source_checkout() -> None:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    commands = [
        (["reproduce_paper.py", "-q"], "All registered manuscript numerics reproduce."),
        (["verify_paper_numerics.py"], "29 confirmed, 0 slip, 0 unknown"),
    ]
    for arguments, expected in commands:
        proc = subprocess.run(
            [sys.executable, str(REPO_ROOT / "scripts" / arguments[0]), *arguments[1:]],
            cwd=REPO_ROOT,
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        assert proc.returncode == 0, proc.stdout + proc.stderr
        assert expected in proc.stdout


def test_release_maintenance_scripts_have_side_effect_free_help() -> None:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    manifest_before = MANIFEST.read_bytes()
    scripts = [
        "build_proof_manifest.py",
        "reproduce_paper.py",
        "verify_paper_numerics.py",
        "verify_release.py",
    ]
    for script in scripts:
        proc = subprocess.run(
            [sys.executable, str(REPO_ROOT / "scripts" / script), "--help"],
            cwd=REPO_ROOT,
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        assert proc.returncode == 0, f"{script}: {proc.stderr}"
        assert "usage:" in proc.stdout.lower(), f"{script}: {proc.stdout}"
        assert "All registered manuscript numerics reproduce." not in proc.stdout
        assert MANIFEST.read_bytes() == manifest_before, f"{script} mutated the manifest"


def test_manifest_builder_check_mode_is_nonmutating() -> None:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    manifest_before = MANIFEST.read_bytes()
    proc = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "build_proof_manifest.py"), "--check"],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode in {0, 1}, proc.stdout + proc.stderr
    output = proc.stdout + proc.stderr
    assert "Current:" in output or "Stale:" in output
    assert MANIFEST.read_bytes() == manifest_before


def test_release_environment_enforces_supported_python_and_exact_packages() -> None:
    scripts = REPO_ROOT / "scripts"
    sys.path.insert(0, str(scripts))
    try:
        import verify_release
    finally:
        sys.path.remove(str(scripts))

    report, _ = verify_release.environment_report()
    assert report["python_requirement"] == "^3.11"
    assert report["python_supported"] is True
    assert verify_release._python_requirement_satisfied("^3.11", (3, 11)) is True
    assert verify_release._python_requirement_satisfied("^3.11", (3, 12)) is True
    assert verify_release._python_requirement_satisfied("^3.11", (3, 13)) is True
    assert verify_release._python_requirement_satisfied("^3.11", (3, 10)) is False
    assert verify_release._python_requirement_satisfied("^3.11", (4, 0)) is False
    assert ".python-version" not in {
        item["path"] for item in json.loads(MANIFEST.read_text(encoding="utf-8"))["artifacts"]
    }


def test_release_verifier_build_toolchain_report_is_coherent() -> None:
    scripts = REPO_ROOT / "scripts"
    sys.path.insert(0, str(scripts))
    try:
        import verify_release
    finally:
        sys.path.remove(str(scripts))

    report, ok = verify_release.build_toolchain_report()
    assert ok is True
    assert report["required_poetry"] == "2.3.3"
    assert report["lockfile_poetry"] == "2.3.3"
    assert report["lock_content_hash_matches"] is True
    assert report["lock_content_hash"] == report["expected_lock_content_hash"]
    assert report["poetry_core"] == "2.3.2"
    assert report["build_backend"] == "poetry.core.masonry.api"
