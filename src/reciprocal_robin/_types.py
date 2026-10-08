"""Typing boundary for dynamically typed mpmath values.

``mpmath`` does not provide complete typed interfaces. The package confines that
dynamic boundary to :data:`Scalar` and :data:`NumericInput`; structural values
remain statically typed.
"""

from typing import Any, TypeAlias

Scalar: TypeAlias = Any
NumericInput: TypeAlias = Any
