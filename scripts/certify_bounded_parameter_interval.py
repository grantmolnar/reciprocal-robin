#!/usr/bin/env python3
"""Generate a bounded-parameter certificate from a source checkout."""

from __future__ import annotations

from release_support import DATA_DIR  # adds ``src/`` to sys.path
from reciprocal_robin.certification.bounded_generate import main


if __name__ == "__main__":
    raise SystemExit(main(default_output_dir=DATA_DIR, default_kappa="1000"))
