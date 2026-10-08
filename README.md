# reciprocal-robin

`reciprocal-robin` contains the numerical code and proof certificates for the
companion paper *Reciprocal LCM powers and Robin's criterion*.
The paper is distributed separately through arXiv.

Repository: <https://github.com/grantmolnar/reciprocal-robin>

The repository has three parts. The ordinary package evaluates reciprocal LCM
powers and their normalized Robin ratios. The exploratory routines study
unnormalized reciprocal-CA candidates. The certification code proves the finite
claims used in the paper by directed interval arithmetic. Its generators and
verifiers are kept separate so that the verifier does not inherit the
calculation it is meant to check.

## Install the numerical package

With Poetry:

```bash
poetry sync
poetry run reciprocal-robin pointwise 5040 2
```

Or from a source checkout with `pip`:

```bash
python -m pip install .
reciprocal-robin pointwise 5040 2
```

Run `reciprocal-robin --help` for the command-line interface. Proof verification
uses the additional pinned tools described below.

## Verify the proof artifacts

Verification is run from a source checkout. The wheel contains the numerical
package, but not the certificates, release scripts, or tests.

The verification dependencies are pinned in `poetry.lock` and
`requirements-verification.txt`. This release supports Python `^3.11`.
Release builds use Poetry 2.3.3 with `poetry-core==2.3.2`.

With Poetry:

```bash
python --version
poetry sync
poetry run python scripts/verify_release.py --strict-environment
```

With a standard virtual environment:

```bash
python -m venv .venv
# POSIX: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements-verification.txt
python scripts/verify_release.py --strict-environment
```

Strict verification requires a supported Python interpreter and the exact
package versions in `requirements-verification.txt`. It checks the certificate and release-file
hashes, Ruff formatting and lint, static typing, both compact certificates, the
numerical values registered from the paper, and the test suite. The exact
manifest hashes belong to release verification rather than ordinary pytest. If
a source edit changes those hashes, rebuild the manifest with

```bash
python scripts/build_proof_manifest.py
```

before rerunning the release verifier. The approximately 3.38-GiB frontier
certificate is regenerated only on explicit request.

## Proof manifest and data

`PROOF_MANIFEST.json` is the byte-exact release inventory: it binds the active
source, verification scripts, environment pins, release documentation, and
mathematical data to the published software release. The `data/` directory has
a narrower role: it contains checked-in certificates and compact certification
records. Its main files are:

- `data/kappa_five_interval_certificate.json`: the finite certificate used for
  the `1 < kappa <= 5` argument and its explicit estimates;
- `data/CERTIFICATION_RESULTS.json`: the compact record for the advertised
  `kappa = 100000` computation; and
- `data/kappa_1000_endpoint_certificate.json`: a small regression certificate
  for the generic bounded-parameter generator and verifier.

The full frontier certificate is too large for the repository. Its SHA-256
digest and certified summary are recorded in `data/CERTIFICATION_RESULTS.json`.
The repository contains the deterministic generator and the independent
verifier needed to reproduce and check it.

See [`CERTIFICATION.md`](CERTIFICATION.md) for the proof boundary, regeneration
commands, and larger-parameter sweeps.

## Package scope

The supported package contains reciprocal LCM-power arithmetic, normalized
Robin ratios, numerical fields and precision controls, normalized extremizers,
and the reusable certification implementations.

The reciprocal-CA routines under `reciprocal_robin.exploration` are exploratory.
They traverse unnormalized extremal candidates. They are not an all-integer
certificate for the main theorem.

Importing `reciprocal_robin` does not change global `mpmath` precision. Code that
needs a particular field or precision should select it explicitly:

```python
from reciprocal_robin import MpfField, using_field, working_precision

with working_precision(50), using_field(MpfField()):
    ...
```

See [`ARCHITECTURE.md`](ARCHITECTURE.md) for the implementation map and the
verifier dependency rules.

## Build a release

The build toolchain is pinned to Poetry 2.3.3 and `poetry-core==2.3.2`. From a
clean checkout, run

```bash
poetry build
```

The wheel contains the runtime package. The sdist also contains the certificates,
verification scripts, and tests. A successful build does not verify the proof;
run `verify_release.py --strict-environment` separately.

## Development

The main checks are:

```bash
python -m mypy src scripts
python -m ruff format --check src scripts tests
python -m ruff check src scripts tests
python scripts/reproduce_paper.py -q
python -m pytest -q -m "not stress"
python -m pytest -q -m stress
```

## Paper and citation

The companion paper is distributed separately through arXiv.

If you use the software, please cite the software release and, when relevant,
the companion paper. `CITATION.cff` contains the software citation metadata.

## License

The software and certification materials are released under the
[MIT License](LICENSE).
