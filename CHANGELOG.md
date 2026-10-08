# Changelog

## 1.0.0 — initial public release

### Certification

- Added a reusable fixed-parameter normalized-endpoint certificate generator and
  an independent verifier.
- Added progressive and explicit-value sweeps for extending the verified
  parameter range on larger machines.
- Included the compact certificate for `1 < kappa <= 5` and the `kappa=1000`
  regression certificate. The larger advertised frontier is represented by its
  certified run record, exact sweep ledger, and certificate digest.
- Scoped CPython's large-integer decimal-conversion override to certification
  operations.
- Kept bounded-generator decisions on live directed intervals rather than on
  rounded decimal diagnostics.
- Strengthened both independent verifiers to reject inconsistent certificate
  metadata before accepting the reconstructed calculation.

### Package and verification

- Separated exploratory reciprocal-CA routines from the supported package root.
- Kept package import neutral with respect to `mpmath` precision and field state.
- Centralized JSON serialization and public argument validation.
- Added release verification for file hashes, exact dependency pins, static
  analysis, certificates, manuscript numerics, and tests.
- Tested the generator/verifier dependency boundary explicitly.
- Removed unreleased compatibility aliases before the first public API freeze.
- Made sampled reciprocal-CA reports respect the requested parsing precision,
  retain both interval bounds, and reject empty or vacuous grids.
- Made certification command dispatch lazy, so invoking a verifier does not
  import its generator.
- Required a successful independent verification before a bounded-parameter
  sweep records telemetry or advances.
- Bound both the JSONL frontier ledger and its CSV mirror to the release record.
- Added a nonmutating manifest check and side-effect-free `--help` behavior to
  the release-maintenance scripts.
- Pinned the release build tools to Poetry 2.3.3 and `poetry-core==2.3.2` and
  check the Poetry lock hash during release verification.
- Separated ordinary tests from byte-exact release integrity: pytest checks
  structure and behavior, while `verify_release.py` checks release hashes.
- Support independent verification on Python `^3.11` with exact pinned
  verification dependencies.

### Distribution

- Split the public software and certification repository from the separately
  distributed arXiv manuscript sources.
- Added MIT licensing and software citation metadata.
