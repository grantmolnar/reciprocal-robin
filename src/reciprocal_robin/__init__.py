"""Numerical tools for reciprocal LCM-power Robin problems.

Importing the package leaves mpmath precision unchanged. Use
:func:`working_precision` to choose precision and :func:`using_field` to choose
point or interval arithmetic. Reciprocal-CA sampling tools live in
:mod:`reciprocal_robin.exploration`.
"""

from importlib.metadata import PackageNotFoundError, version

from ._field import MpfField, MpiField, RealField, Trilean, current_field, using_field
from ._precision import DEFAULT_DPS, set_default_precision, working_precision
from .arithmetic import F_local, rksigma_pp, sigma_minus_one_pp
from .factorization import Factorization
from .normalized import (
    boundary_difference,
    exceptional_set,
    kappa_star,
    rho_hat_one,
    rho_k,
    robin_ratio,
)
from .normalized_ca import normalized_local_objective, rho_pp, square_entry_gap_k5

try:
    __version__ = version("reciprocal-robin")
except PackageNotFoundError:  # uninstalled source checkout
    __version__ = "0+unknown"

__all__ = [
    "DEFAULT_DPS",
    "F_local",
    "Factorization",
    "MpfField",
    "MpiField",
    "RealField",
    "Trilean",
    "boundary_difference",
    "current_field",
    "exceptional_set",
    "kappa_star",
    "normalized_local_objective",
    "rho_hat_one",
    "rho_k",
    "rho_pp",
    "rksigma_pp",
    "robin_ratio",
    "set_default_precision",
    "sigma_minus_one_pp",
    "square_entry_gap_k5",
    "using_field",
    "working_precision",
]
