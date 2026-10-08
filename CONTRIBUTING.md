# Contributing

This repository contains both ordinary numerical software and code used in a
computer-assisted proof. Ordinary refactors should satisfy the usual software
checks. A change on the proof path must also preserve the meaning and
independence of the certificates.

## Development setup

Use Poetry 2.3.3, as required by `pyproject.toml` and recorded in the lockfile.
The PEP 517 backend is pinned separately to `poetry-core==2.3.2`.

```bash
poetry --version
python --version  # any Python satisfying ^3.11
poetry sync
```

Run:

```bash
poetry run python -m mypy src scripts
poetry run python -m ruff format --check src scripts tests
poetry run python -m ruff check src scripts tests
poetry run python scripts/reproduce_paper.py -q
poetry run python -m pytest -q -m "not stress"
poetry run python -m pytest -q -m stress
```

For a release candidate, format first, then rebuild the byte-exact manifest,
and finally run the release verifier:

```bash
poetry run python -m ruff format src scripts tests
poetry run python -m ruff check --fix src scripts tests
poetry run python scripts/build_proof_manifest.py
poetry run python scripts/verify_release.py --strict-environment
```

Formatting changes the manifest hashes. This is a release-integrity issue, not a
behavioral-test failure: pytest checks the manifest's structure and inventory,
while `verify_release.py` checks its exact hashes.

## Test organization

The test tree is organized by the kind of assurance each check provides:

- `tests/unit/` covers deterministic API, arithmetic, state, CLI, and walker behavior.
- `tests/properties/` contains Hypothesis tests for structural and theorem-level invariants.
- `tests/certification/` exercises certificate generation, independent verification, and
  published certification records.
- `tests/release/` checks packaging, manifest, environment, and source-checkout contracts.
- `tests/paper/` ties registered manuscript numerics to the code.

The `stress` marker is orthogonal to this layout: slower tests remain beside the subsystem
they exercise and are selected with `-m stress`. Repository-relative tests should use the
shared `tests._repo.REPO_ROOT` helper rather than infer the checkout root from their own
directory depth.

## Changes on the proof path

Do not weaken generator/verifier independence. In particular, the bounded
verifier must not import the bounded generator, `reciprocal_robin.arithmetic`,
or `reciprocal_robin.normalized_ca`. The `kappa <= 5` verifier keeps its own
theorem-level interval formulas for the same reason.

A change to a mathematical formula, endpoint or event definition, certificate
field meaning, or analytic handoff is a proof change. It requires a fresh
technical review. A performance change should leave certificate semantics
unchanged and should be checked by deterministic regeneration followed by the
independent verifier.

A larger verified parameter does not require a source change. Use the generic
certificate or sweep commands in `CERTIFICATION.md`, review the resulting
artifact, and change the advertised frontier only after independent
verification.

## Public API

Keep new public interfaces small, typed, and documented. Add them when they
serve a concrete use of the package. Exploratory routines belong under
`reciprocal_robin.exploration` unless they become part of the supported
numerical or certification API.

This is the first public release. There are no public legacy clients to support,
so obsolete private-development names should be removed rather than preserved
through compatibility aliases.

## Generated files

Do not commit ad hoc sweep directories, large generated frontier certificates,
Python caches, or local logs. Compact certificates and release metadata are updated through the repository
scripts.

## Release packaging

Build both distributions from a clean checkout with the pinned toolchain:

```bash
poetry build
```

Inspect both archives before upload. The wheel should contain the importable
package but not repository-only certificate data or tests. The sdist should
contain the proof-verification material declared in `pyproject.toml`.

Packaging and proof verification are separate checks. The proof is verified
from a source checkout with `verify_release.py --strict-environment`.
