"""Verify a release from a clean checkout."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import re
import subprocess
import sys
import tomllib
from typing import Any

from release_support import (
    PROOF_MANIFEST,
    REPO_ROOT,
    REQUIREMENTS_FILE,
    read_json,
    sha256_file,
)
from reciprocal_robin._serialization import json_text

MANIFEST_PATH = PROOF_MANIFEST
REQUIREMENTS_PATH = REQUIREMENTS_FILE


def parse_pins() -> dict[str, str]:
    """Parse the exact release pins, rejecting malformed or duplicate entries."""
    pins: dict[str, str] = {}
    for line_number, raw in enumerate(
        REQUIREMENTS_PATH.read_text(encoding="utf-8").splitlines(), start=1
    ):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        match = re.fullmatch(r"([A-Za-z0-9_.-]+)==([^\s]+)", line)
        if match is None:
            raise ValueError(
                f"invalid exact pin on line {line_number} of {REQUIREMENTS_PATH.name}: {line!r}"
            )
        name = match.group(1).lower().replace("_", "-")
        if name in pins:
            raise ValueError(f"duplicate verification pin for {name!r}")
        pins[name] = match.group(2)
    if not pins:
        raise ValueError(f"{REQUIREMENTS_PATH.name} contains no verification pins")
    return pins


def _python_requirement() -> str:
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        requirement = tomllib.load(handle)["tool"]["poetry"]["dependencies"]["python"]
    if not isinstance(requirement, str):
        raise ValueError("tool.poetry.dependencies.python must be a string")
    return requirement


def _python_requirement_satisfied(requirement: str, version: tuple[int, int] | None = None) -> bool:
    # This parser implements the current ^3.11 requirement exactly. A changed
    # project constraint must be handled here rather than interpreted implicitly.
    if requirement != "^3.11":
        raise ValueError(f"unsupported Python requirement syntax: {requirement!r}")
    major_minor = (sys.version_info.major, sys.version_info.minor) if version is None else version
    return (3, 11) <= major_minor < (4, 0)


def environment_report() -> tuple[dict[str, Any], bool]:
    rows = []
    packages_exact = True
    for name, expected in parse_pins().items():
        try:
            installed = importlib.metadata.version(name)
            status = "exact" if installed == expected else "drift"
        except importlib.metadata.PackageNotFoundError:
            installed = None
            status = "missing"
        packages_exact &= status == "exact"
        rows.append(
            {
                "package": name,
                "expected": expected,
                "installed": installed,
                "status": status,
            }
        )

    installed_python = sys.version.split()[0]
    python_requirement = _python_requirement()
    python_supported = _python_requirement_satisfied(python_requirement)
    strict_ok = packages_exact and python_supported
    return {
        "python": installed_python,
        "python_requirement": python_requirement,
        "python_supported": python_supported,
        "packages_exact": packages_exact,
        "strict_environment_ok": strict_ok,
        "packages": rows,
    }, strict_ok


def _citation_version() -> str:
    for line in (REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8").splitlines():
        if line.startswith("version:"):
            return line.split(":", 1)[1].strip()
    raise ValueError("CITATION.cff has no top-level version")


def _pyproject_version() -> str:
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        return str(tomllib.load(handle)["tool"]["poetry"]["version"])


def _poetry_233_content_hash(pyproject: dict[str, Any]) -> str:
    """Reproduce Poetry 2.3.3's lockfile content hash."""
    legacy_keys = ["dependencies", "source", "extras", "dev-dependencies"]
    relevant_keys = [*legacy_keys, "group"]
    relevant_project_keys = ["requires-python", "dependencies", "optional-dependencies"]

    project_content = pyproject.get("project", {})
    group_content = pyproject.get("dependency-groups", {})
    tool_poetry_content = pyproject.get("tool", {}).get("poetry", {})

    relevant_project_content = {
        key: project_content[key] for key in relevant_project_keys if key in project_content
    }
    relevant_poetry_content: dict[str, Any] = {}
    for key in relevant_keys:
        data = tool_poetry_content.get(key)
        if data is None and (key not in legacy_keys or relevant_project_content or group_content):
            continue
        relevant_poetry_content[key] = data

    relevant_content: dict[str, Any]
    if relevant_project_content or group_content:
        relevant_content = {}
        if relevant_project_content:
            relevant_content["project"] = relevant_project_content
        if group_content:
            relevant_content["dependency-groups"] = group_content
        relevant_content["tool"] = {"poetry": relevant_poetry_content}
    else:
        relevant_content = relevant_poetry_content

    encoded = json.dumps(relevant_content, sort_keys=True).encode()
    return hashlib.sha256(encoded).hexdigest()


def build_toolchain_report() -> tuple[dict[str, Any], bool]:
    """Check the pinned build tools against the lockfile."""
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        pyproject = tomllib.load(handle)

    poetry_requirement = pyproject.get("tool", {}).get("poetry", {}).get("requires-poetry")
    poetry_match = (
        re.fullmatch(r"==([^\s]+)", poetry_requirement)
        if isinstance(poetry_requirement, str)
        else None
    )
    required_poetry = poetry_match.group(1) if poetry_match else None

    lock_header = (REPO_ROOT / "poetry.lock").read_text(encoding="utf-8").splitlines()[0]
    lock_pattern = (
        r"# This file is automatically @generated by Poetry ([^ ]+) "
        r"and should not be changed by hand\."
    )
    lock_match = re.fullmatch(lock_pattern, lock_header)
    lock_poetry = lock_match.group(1) if lock_match else None
    with (REPO_ROOT / "poetry.lock").open("rb") as handle:
        lock_data = tomllib.load(handle)
    lock_content_hash = lock_data.get("metadata", {}).get("content-hash")
    expected_content_hash = _poetry_233_content_hash(pyproject)
    lock_fresh = lock_content_hash == expected_content_hash

    build_system = pyproject.get("build-system", {})
    build_requires = build_system.get("requires")
    core_version = None
    if isinstance(build_requires, list) and len(build_requires) == 1:
        requirement = build_requires[0]
        if isinstance(requirement, str):
            core_match = re.fullmatch(r"poetry-core==([^\s]+)", requirement)
            if core_match:
                core_version = core_match.group(1)

    backend = build_system.get("build-backend")
    ok = (
        required_poetry is not None
        and required_poetry == lock_poetry
        and lock_fresh
        and core_version is not None
        and backend == "poetry.core.masonry.api"
    )
    return {
        "required_poetry": required_poetry,
        "lockfile_poetry": lock_poetry,
        "poetry_matches_lockfile": required_poetry is not None and required_poetry == lock_poetry,
        "lock_content_hash": lock_content_hash,
        "expected_lock_content_hash": expected_content_hash,
        "lock_content_hash_matches": lock_fresh,
        "poetry_core": core_version,
        "build_backend": backend,
        "exact": ok,
    }, ok


def verify_manifest() -> tuple[dict[str, Any], bool]:
    """Check manifest hashes and reconstruct the expected inventory."""
    from build_proof_manifest import build_manifest

    manifest = read_json(MANIFEST_PATH)
    if not isinstance(manifest, dict):
        raise ValueError("proof manifest must be a JSON object")
    expected = build_manifest()
    inventory_matches = manifest == expected

    rows = []
    metadata_versions = {
        "manifest": manifest.get("package_version"),
        "pyproject": _pyproject_version(),
        "citation": _citation_version(),
    }
    metadata_matches = len(set(metadata_versions.values())) == 1
    ok = metadata_matches and inventory_matches
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list):
        raise ValueError("proof manifest artifacts must be a list")
    for item in artifacts:
        if not isinstance(item, dict) or not isinstance(item.get("path"), str):
            raise ValueError("every proof manifest artifact must be an object with a path")
        path = REPO_ROOT / item["path"]
        exists = path.is_file()
        actual = sha256_file(path) if exists else None
        matches = exists and actual == item.get("sha256")
        ok &= matches
        rows.append({"path": item["path"], "exists": exists, "hash_matches": matches})
    return {
        "schema_version": manifest.get("schema_version"),
        "mathematical_release": manifest.get("mathematical_release"),
        "package_versions": metadata_versions,
        "package_versions_match": metadata_matches,
        "inventory_matches_builder": inventory_matches,
        "artifacts": rows,
    }, ok


def run_step(label: str, args: list[str], env: dict[str, str]) -> bool:
    print(f"\n=== {label} ===", flush=True)
    proc = subprocess.run(args, cwd=REPO_ROOT, env=env)
    if proc.returncode:
        print(f"FAILED: {label} (exit {proc.returncode})", file=sys.stderr)
        return False
    print(f"PASSED: {label}")
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--strict-environment", action="store_true")
    parser.add_argument("--skip-tests", action="store_true")
    args = parser.parse_args()
    if sys.version_info < (3, 11):
        print("Python 3.11 or newer is required.", file=sys.stderr)
        return 2

    try:
        env_data, env_ok = environment_report()
    except ValueError as exc:
        print(f"Invalid verification environment specification: {exc}", file=sys.stderr)
        return 2
    print("=== Environment ===")
    print(json_text(env_data))
    if args.strict_environment and not env_ok:
        print(
            "Environment does not satisfy the supported Python range and exact dependency pins.",
            file=sys.stderr,
        )
        return 2

    try:
        build_data, build_ok = build_toolchain_report()
    except (IndexError, KeyError, TypeError, ValueError) as exc:
        print(f"Invalid release build-toolchain metadata: {exc}", file=sys.stderr)
        return 2
    print("\n=== Release build toolchain ===")
    print(json_text(build_data))
    if not build_ok:
        print(
            "Release build-toolchain metadata is not exact or internally consistent.",
            file=sys.stderr,
        )
        return 1

    try:
        manifest_data, manifest_ok = verify_manifest()
    except (KeyError, TypeError, ValueError) as exc:
        print(f"Invalid proof manifest: {exc}", file=sys.stderr)
        return 2
    print("\n=== Proof-artifact hashes ===")
    print(json_text(manifest_data))
    if not manifest_ok:
        print("Proof-artifact manifest verification failed.", file=sys.stderr)
        return 1

    env = os.environ.copy()
    src = str(REPO_ROOT / "src")
    env["PYTHONPATH"] = src if not env.get("PYTHONPATH") else src + os.pathsep + env["PYTHONPATH"]
    steps = [
        (
            "strict static typing",
            [sys.executable, "-m", "mypy", "src", "scripts"],
        ),
        (
            "Ruff formatting",
            [sys.executable, "-m", "ruff", "format", "--check", "src", "scripts", "tests"],
        ),
        (
            "Ruff lint",
            [sys.executable, "-m", "ruff", "check", "src", "scripts", "tests"],
        ),
        (
            "kappa<=5 finite theorem certificate",
            [sys.executable, "scripts/verify_kappa_five_certificate.py"],
        ),
        (
            "kappa=1000 bounded-certificate regression fixture",
            [sys.executable, "scripts/verify_bounded_parameter_certificate.py"],
        ),
        ("registered manuscript numerics", [sys.executable, "scripts/reproduce_paper.py", "-q"]),
    ]
    if not args.skip_tests:
        steps.extend(
            [
                (
                    "active test suite",
                    [sys.executable, "-m", "pytest", "-q", "-m", "not stress"],
                ),
                (
                    "stress test suite",
                    [sys.executable, "-m", "pytest", "-q", "-m", "stress"],
                ),
            ]
        )

    all_ok = True
    for label, command in steps:
        if not run_step(label, command, env):
            all_ok = False
            break
    print("\n=== Release verification summary ===")
    print(f"artifact_hashes: {'PASS' if manifest_ok else 'FAIL'}")
    print(f"build_toolchain: {'PASS' if build_ok else 'FAIL'}")
    print(f"environment: {'PASS' if env_ok else 'DRIFT (non-fatal)'}")
    print(f"proof_and_tests: {'PASS' if all_ok else 'FAIL'}")
    return 0 if manifest_ok and all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
