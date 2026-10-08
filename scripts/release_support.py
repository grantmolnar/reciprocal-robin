"""Shared repository paths and serialization for certification and release scripts."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from reciprocal_robin._serialization import read_json, write_json  # noqa: E402

DATA_DIR = REPO_ROOT / "data"
PROOF_MANIFEST = REPO_ROOT / "PROOF_MANIFEST.json"
REQUIREMENTS_FILE = REPO_ROOT / "requirements-verification.txt"

KAPPA_FIVE_CERTIFICATE = DATA_DIR / "kappa_five_interval_certificate.json"
KAPPA_1000_CERTIFICATE = DATA_DIR / "kappa_1000_endpoint_certificate.json"
CERTIFICATION_RESULTS = DATA_DIR / "CERTIFICATION_RESULTS.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = [
    "REPO_ROOT",
    "SRC_DIR",
    "DATA_DIR",
    "PROOF_MANIFEST",
    "REQUIREMENTS_FILE",
    "KAPPA_FIVE_CERTIFICATE",
    "KAPPA_1000_CERTIFICATE",
    "CERTIFICATION_RESULTS",
    "read_json",
    "write_json",
    "sha256_file",
]
