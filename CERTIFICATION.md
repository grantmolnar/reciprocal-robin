# Certification

The certification code supplies the finite computer-assisted part of the
bounded-parameter results in *Reciprocal LCM powers and Robin's criterion*. The analytic reduction is proved in the paper; the code certifies
the finite inequalities left by that reduction.

## Independent verification

There are two certificate families, each with a generator and a separate
verifier.

- `reciprocal_robin.certification.kappa_five_generate` and
  `kappa_five_verify` handle the finite certificate for `1 < kappa <= 5` and
  the explicit Mertens inputs.
- `reciprocal_robin.certification.bounded_generate` and `bounded_verify` handle
  the generic fixed-parameter normalized-endpoint certificate.

The verifiers do not import their generators or the package's theorem-specific
arithmetic. The bounded verifier reconstructs the finite tail, local hulls,
endpoint transitions, and Robin margins directly from the certificate
witnesses. The `kappa <= 5` verifier reimplements the prime-power formulas and
reconstructs the complete compact artifact. Generic infrastructure such as
serialization, interval arithmetic, prime enumeration, and command-line
orchestration may be shared.

This separation is part of the proof. A performance refactor must not replace
independent recomputation by a call into the generator.

The bounded verifier also checks redundant discrete data against its own
reconstruction: local hull order, global transition order, event metadata,
endpoint values, and the two `all_positive` flags. For the compact
`kappa <= 5` certificate, the independently reconstructed JSON object must
agree with the checked-in object after all finite checks pass.

### Decimal diagnostics

Historical fields ending in `_lo`, `_hi`, or `_lower` are decimal renderings of
binary interval endpoints. They are useful for inspection and run-to-run
comparison, but the bounded 1.0 and `kappa <= 5` 1.1 schemas write them with
`mpmath.nstr`. They are therefore rounded decimal strings, not independently
outward-rounded decimal enclosures.

No proof decision depends on those strings. Generators make sign and ordering
decisions with the directed interval endpoints before serialization, and the
verifiers recompute the intervals from the discrete witnesses and manuscript
formulas. The bounded generator likewise retains the entry and event intervals
through ordering and separation rather than reparsing their decimal renderings.

Changing the representation would change the digest of the existing large
frontier certificate. An outward-rounded textual representation should
therefore be introduced, if wanted, as a new certificate schema.

The `symbolic_square_entry_identity` field in the compact `kappa <= 5`
certificate records the exact algebraic square-entry lemma from the manuscript.
It is provenance for that symbolic argument, not an additional numerical
interval check.

## Released artifacts

`PROOF_MANIFEST.json` is the byte-exact release inventory. It records the role,
verifier, size, and SHA-256 digest of the proof data, environment pins,
verification documentation, release scripts, and active package source. The
`data/` directory itself is reserved for certificates and compact certification
records. The release verifier reconstructs the expected
inventory from the checkout before comparing hashes, so the manifest cannot
silently omit a newly added package module or verification script.

The current release records `100000` as the advertised parameter in
`data/CERTIFICATION_RESULTS.json`. Its full endpoint certificate is
3,624,666,835 bytes and is not stored in Git. The compact record gives its
SHA-256 digest, verification environment, finite-tail handoff, endpoint count,
and certified margins. The certificate is reproduced from the released
generator and checked by the independent verifier.

The checked-in `kappa_1000_endpoint_certificate.json` runs through the same
generic generator and verifier at a size suitable for routine regression tests.

## Verify the release

After installing the pinned environment, run:

```bash
python scripts/verify_release.py --strict-environment
```

To check the artifacts without running the test suite:

```bash
python scripts/verify_release.py --strict-environment --skip-tests
```

Strict mode requires an interpreter in the Python range declared by
`pyproject.toml` and exact matches for every package in
`requirements-verification.txt`.

## Regenerate the compact certificates

For the `kappa <= 5` certificate:

```bash
python scripts/certify_kappa_five_interval.py
python scripts/verify_kappa_five_certificate.py
```

For the compact generic regression certificate:

```bash
python scripts/certify_bounded_parameter_interval.py --kappa 1000
python scripts/verify_bounded_parameter_certificate.py
```

## Reproduce the advertised frontier

Generate the large certificate outside the repository tree and verify it
independently:

```bash
python scripts/certify_bounded_parameter_interval.py \
    --kappa 100000 \
    --output kappa_100000_endpoint_certificate.json
python scripts/verify_bounded_parameter_certificate.py \
    --certificate kappa_100000_endpoint_certificate.json
```

For the current release its SHA-256 digest is
`d94f3499a46c155fbcdd7e0a6f20f773ca8dd031f2af9b4217246aff9f5456c3`.
Timing depends on the machine; the digest and certified inequalities do not.

## Extend the parameter range

No source change is needed to try a larger integer parameter. For one value:

```bash
python scripts/certify_bounded_parameter_interval.py \
    --kappa 2000 \
    --output kappa_2000.json
python scripts/verify_bounded_parameter_certificate.py \
    --certificate kappa_2000.json
```

To search a sequence automatically:

```bash
python scripts/sweep_bounded_parameter.py \
    --start 2000 \
    --factor 2 \
    --timeout-minutes 30
```

Resume an interrupted sweep with `--resume`, or give selected values with
`--values`:

```bash
python scripts/sweep_bounded_parameter.py \
    --values 120000,140000,160000 \
    --timeout-minutes 960 \
    --output-dir sweep-refine
```

The sweep advances only after the independent verifier succeeds. A timeout,
out-of-memory failure, or generator error marks a computational frontier, not a
mathematical counterexample.

If an interval comparison is indeterminate, retrying at higher precision may
separate its sign. Higher precision cannot turn a genuinely negative margin
into a proof.

The bounded generator and verifier temporarily relax CPython's decimal
integer-conversion guard while processing exact endpoint witnesses, then
restore the previous setting. This guard is process-global, not thread-local.
Certification should therefore not run concurrently in one interpreter with
code that relies on the default denial-of-service limit. No global
`PYTHONINTMAXSTRDIGITS` setting is needed.

## Installed-package commands

The same operations are available from an installed package:

```bash
reciprocal-robin-certify-bounded --kappa 2000 --output kappa_2000.json
reciprocal-robin-verify-bounded --certificate kappa_2000.json
reciprocal-robin-sweep-bounded --start 2000 --factor 2
```

The module dispatcher provides shorter command names:

```bash
python -m reciprocal_robin.certification generate --kappa 2000 --output kappa_2000.json
python -m reciprocal_robin.certification verify --certificate kappa_2000.json
python -m reciprocal_robin.certification sweep --start 2000 --factor 2
python -m reciprocal_robin.certification generate-kappa-five --output kappa_five.json
python -m reciprocal_robin.certification verify-kappa-five --certificate kappa_five.json
```

## Promoting a stronger frontier

A completed sweep is not by itself enough to change the advertised frontier.
Before changing the paper or release metadata:

1. retain the exact sweep ledger and certificate digest;
2. rerun the certificate deterministically when practical;
3. verify it independently in the pinned environment;
4. inspect the finite-tail handoff and certified margins; and
5. update the compact run record and proof manifest.

After changing a proof artifact or release file, rebuild the manifest and check
the release:

```bash
python scripts/build_proof_manifest.py
python scripts/verify_release.py --strict-environment
```
