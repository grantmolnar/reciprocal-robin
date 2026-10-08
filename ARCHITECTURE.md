# Architecture

There are three layers: ordinary numerical code, exploratory extremal code, and
proof certification. The third is separated because an independent verifier
should not inherit theorem-specific calculations from the generator it checks.

## Numerical foundation

`_field.py` provides point (`MpfField`) and interval (`MpiField`) arithmetic.
`_precision.py` scopes changes to numerical precision, so importing the package
does not change global `mpmath` state. `_validation.py` checks mathematical
domains at public interfaces. `_serialization.py` provides atomic UTF-8 JSON
I/O.

`arithmetic.py`, `factorization.py`, and `normalized.py` implement the
prime-power formulas, their multiplicative assembly, normalized reciprocal LCM
powers, Robin ratios, and finite exceptional-set calculations.

## Extremal calculations

`normalized_ca.py` implements the local extremizers used in the
bounded-parameter certificate.

The unnormalized reciprocal-CA calculations live in `ca.py`, `successor.py`,
`walker.py`, `incremental_ratio.py`, and `ca_grid.py`. Their supported public
interface is `reciprocal_robin.exploration`. These routines are useful for
pointwise exploration, but they are not used as an all-integer certificate for
the main theorem.

## Certification

`reciprocal_robin.certification` contains two certificate families:

- `kappa_five_generate.py` and `kappa_five_verify.py` certify the finite inputs
  to the `1 < kappa <= 5` argument;
- `bounded_generate.py` and `bounded_verify.py` handle the generic
  fixed-parameter normalized-endpoint argument.

`sweep.py` only coordinates runs. It starts the generator and verifier in fresh
processes, records their outputs, and advances to another parameter only after
the verifier succeeds.

The commands in `scripts/` are thin wrappers around these modules. The unified
`python -m reciprocal_robin.certification` dispatcher imports only the selected
command, so running a verifier does not import its generator as a side effect of
command dispatch.

## Verifier independence

The bounded verifier does not import its generator,
`reciprocal_robin.arithmetic`, or `reciprocal_robin.normalized_ca`. It
reconstructs the theorem-specific formulas independently. The `kappa <= 5`
verifier follows the same rule.

Generic infrastructure may be shared when it carries no theorem-specific
calculation: serialization, interval backends, prime enumeration, process
coordination, and immutable certificate formats.

## Paper numerics

`paper_numerics.py` records the pointwise numerical values printed in the paper. `scripts/reproduce_paper.py` and
`scripts/verify_paper_numerics.py` check them without requiring the TeX sources.

## Release verification

`PROOF_MANIFEST.json` is the release inventory, not another certificate-data
directory. It records the files fixed by the release and their hashes; checked-in
certificates and compact run records remain under `data/`.
`scripts/verify_release.py` first reconstructs the expected inventory from the
checkout, then checks the hashes, compact certificates, paper numerics, static
analysis, and tests. A newly added package or proof script therefore cannot be
silently omitted from an old manifest.

The active Python package is hashed recursively, so the release hashes bind the
certification code as well as the certificate data.

See `CERTIFICATION.md` for the certificate procedures and `CONTRIBUTING.md` for
the rules governing changes to the certified calculations.
