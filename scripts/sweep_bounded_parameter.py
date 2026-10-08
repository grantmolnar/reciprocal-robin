#!/usr/bin/env python3
"""Run a bounded-parameter certification sweep from a source checkout."""

from __future__ import annotations

from release_support import REPO_ROOT  # adds ``src/`` to sys.path
from reciprocal_robin.certification.sweep import main


if __name__ == "__main__":
    raise SystemExit(main(default_repo=REPO_ROOT))
