#!/usr/bin/env python3
"""Generate the ``kappa <= 5`` finite certificate from a source checkout."""

from __future__ import annotations

from release_support import KAPPA_FIVE_CERTIFICATE  # adds ``src/`` to sys.path
from reciprocal_robin.certification.kappa_five_generate import main


if __name__ == "__main__":
    raise SystemExit(main(default_output=KAPPA_FIVE_CERTIFICATE))
