"""Certificate generators and independent verifiers.

``kappa_five_generate`` and ``kappa_five_verify`` handle the finite inputs for
``1 < kappa <= 5``. ``bounded_generate`` and ``bounded_verify`` handle the
generic fixed-parameter normalized-endpoint certificate.

The verifiers do not import their generators or the package's theorem-specific
arithmetic. ``sweep`` repeats the bounded generation-and-verification pair over
a sequence of integer parameters.
"""

__all__: list[str] = []
