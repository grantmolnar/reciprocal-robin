"""Reciprocal-CA construction, successor walks, and sampled-grid calculations.

These routines study the extremal candidates themselves. A sampled CA walk does
not control the integers between consecutive candidates; theorem certificates
therefore use the separate certification package.
"""

from .ca import MaxExponentExceeded, ThresholdTooLarge, reciprocal_ca, x_threshold
from .incremental_ratio import (
    IncrementalRatioState,
    IncrementalRatioWalker,
    IncrementalStep,
    SuccessorInvariantError,
)
from .successor import CASuccessorEnumerator
from .ca_grid import (
    CAWalkReport,
    CertifiedStep,
    audit_ca_grid,
    print_step,
    max_robin_ratio_on_ca_grid,
)
from .walker import CAStep, walk_ca

__all__ = [
    "CAStep",
    "CASuccessorEnumerator",
    "CAWalkReport",
    "CertifiedStep",
    "IncrementalRatioState",
    "IncrementalRatioWalker",
    "IncrementalStep",
    "MaxExponentExceeded",
    "SuccessorInvariantError",
    "ThresholdTooLarge",
    "audit_ca_grid",
    "print_step",
    "reciprocal_ca",
    "max_robin_ratio_on_ca_grid",
    "walk_ca",
    "x_threshold",
]
