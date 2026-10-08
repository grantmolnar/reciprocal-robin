#!/usr/bin/env python3
"""Independently verify a bounded certificate from a source checkout."""

from __future__ import annotations

from release_support import KAPPA_1000_CERTIFICATE  # adds ``src/`` to sys.path
from reciprocal_robin.certification.bounded_verify import main


if __name__ == "__main__":
    raise SystemExit(main(default_certificate=KAPPA_1000_CERTIFICATE))
