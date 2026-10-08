#!/usr/bin/env python3
"""Independently verify the ``kappa <= 5`` certificate from a source checkout."""

from __future__ import annotations

from release_support import KAPPA_FIVE_CERTIFICATE  # adds ``src/`` to sys.path
from reciprocal_robin.certification.kappa_five_verify import main


if __name__ == "__main__":
    raise SystemExit(main(default_certificate=KAPPA_FIVE_CERTIFICATE))
