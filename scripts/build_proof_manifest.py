"""Build the checked-in inventory of certification and release artifacts."""

from __future__ import annotations

import argparse
import sys
import tomllib
from collections.abc import Sequence
from typing import Any

from release_support import (
    CERTIFICATION_RESULTS,
    PROOF_MANIFEST,
    REPO_ROOT,
    read_json,
    sha256_file,
    write_json,
)

OUTPUT = PROOF_MANIFEST
CERTIFICATION_RECORD = read_json(CERTIFICATION_RESULTS)
ADVERTISED_KAPPA = int(CERTIFICATION_RECORD["advertised_kappa"])
BASE_ARTIFACTS = [
    {
        "path": "data/kappa_five_interval_certificate.json",
        "role": (
            "Directed-interval certificate for the unconditional results. "
            "The interval 1 < kappa <= 5 uses the squared-primorial range, "
            "Mertens extension, normalized endpoint, and low-support margin; "
            "the uniform explicit bound uses the primorial threshold and shares "
            "the Mertens extension; the bounded-parameter analytic handoff "
            "also invokes that certified Mertens extension."
        ),
        "verifier": "scripts/verify_kappa_five_certificate.py",
        "load_bearing": True,
    },
    {
        "path": "data/CERTIFICATION_RESULTS.json",
        "role": (
            f"Certified run record for the advertised all-integer interval "
            f"1 < kappa <= {ADVERTISED_KAPPA}. It records the verification environment, "
            "frontier certificate digest and size, finite-tail handoff, endpoint "
            "count, and certified margins. The frontier certificate itself is "
            "reproducible from this release but is not bundled."
        ),
        "verifier": "scripts/verify_release.py",
        "load_bearing": True,
    },
    {
        "path": "data/kappa_1000_endpoint_certificate.json",
        "role": (
            "Compact regression fixture for the generic bounded-parameter "
            "generator and independent verifier; not the proof artifact for the "
            f"advertised kappa={ADVERTISED_KAPPA} frontier."
        ),
        "verifier": "scripts/verify_bounded_parameter_certificate.py",
        "load_bearing": False,
    },
    {
        "path": "pyproject.toml",
        "role": "Package metadata and direct dependency specification.",
        "verifier": "scripts/verify_release.py",
        "load_bearing": False,
    },
    {
        "path": "CITATION.cff",
        "role": "Archive-ready software citation metadata.",
        "verifier": "scripts/verify_release.py",
        "load_bearing": False,
    },
    {
        "path": "poetry.lock",
        "role": "Poetry dependency lock used for the release environment.",
        "verifier": "scripts/verify_release.py",
        "load_bearing": False,
    },
    {
        "path": "requirements-verification.txt",
        "role": "Portable exact package pins for proof verification without Poetry.",
        "verifier": "scripts/verify_release.py --strict-environment",
        "load_bearing": False,
    },
    {
        "path": "README.md",
        "role": "Public entry point for the software and certificate-verification release.",
        "verifier": "scripts/verify_release.py",
        "load_bearing": False,
    },
    {
        "path": "ARCHITECTURE.md",
        "role": "Public description of the certification trust boundary and package architecture.",
        "verifier": "scripts/verify_release.py",
        "load_bearing": False,
    },
    {
        "path": "CERTIFICATION.md",
        "role": "Public proof-artifact verification and regeneration instructions.",
        "verifier": "scripts/verify_release.py",
        "load_bearing": False,
    },
]


def specifications() -> list[dict[str, Any]]:
    specs = list(BASE_ARTIFACTS)
    proof_scripts = sorted((REPO_ROOT / "scripts").glob("*.py"))
    for script in proof_scripts:
        relative = str(script.relative_to(REPO_ROOT))
        specs.append(
            {
                "path": relative,
                "role": "Certificate generation, verification, or release-support script.",
                "verifier": "scripts/verify_release.py",
                "load_bearing": False,
            }
        )
    for path in sorted((REPO_ROOT / "src" / "reciprocal_robin").rglob("*.py")):
        specs.append(
            {
                "path": str(path.relative_to(REPO_ROOT)),
                "role": "Package source used by the certificate-verification scripts.",
                "verifier": "scripts/verify_release.py",
                "load_bearing": False,
            }
        )
    return specs


def package_version() -> str:
    """Read the Python distribution version from the canonical metadata."""
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        return str(tomllib.load(handle)["tool"]["poetry"]["version"])


def build_manifest() -> dict[str, Any]:
    entries = []
    for spec in specifications():
        path = REPO_ROOT / spec["path"]
        if not path.is_file():
            raise FileNotFoundError(path)
        entries.append({**spec, "sha256": sha256_file(path), "bytes": path.stat().st_size})
    return {
        "schema_version": "proof-artifact-inventory-1.1.0",
        "mathematical_release": "reciprocal-robin-1.0.0",
        "package_version": package_version(),
        "manuscript": "Reciprocal LCM powers and Robin's criterion",
        "verification_command": "python scripts/verify_release.py",
        "strict_environment_command": "python scripts/verify_release.py --strict-environment",
        "artifacts": entries,
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build or check the certificate and release manifest.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="check that the checked-in manifest is current without rewriting it",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    expected = build_manifest()
    if args.check:
        current = read_json(OUTPUT)
        if current == expected:
            print(f"Current: {OUTPUT.relative_to(REPO_ROOT)}")
            return 0
        print(f"Stale: {OUTPUT.relative_to(REPO_ROOT)}", file=sys.stderr)
        return 1

    write_json(OUTPUT, expected)
    print(f"Wrote {OUTPUT.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
