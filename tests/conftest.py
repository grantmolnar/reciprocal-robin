"""Shared deterministic numerical context for the test suite."""

from __future__ import annotations

import pytest

from reciprocal_robin import DEFAULT_DPS, working_precision


@pytest.fixture(autouse=True)
def _release_precision():
    """Run each test at the release baseline without leaking global state."""
    with working_precision(DEFAULT_DPS):
        yield
