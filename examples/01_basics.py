"""A short tour of the reciprocal-robin numerical API.

Run from the repository root with:

    poetry run python examples/01_basics.py
"""

from mpmath import mpf

from reciprocal_robin import (
    Factorization,
    MpfField,
    exceptional_set,
    kappa_star,
    rho_k,
    robin_ratio,
    using_field,
    working_precision,
)
from reciprocal_robin.exploration import reciprocal_ca


def section(title: str) -> None:
    print(f"\n{'=' * 60}\n  {title}\n{'=' * 60}")


def main() -> None:
    with working_precision(50), using_field(MpfField()):
        section("1. Pointwise evaluation at (n, kappa) = (36, 2)")
        print(f"  rho^[2](36)                = {float(rho_k(36, mpf(2))):.10f}")
        print(f"  rho^[2](36) / Robin scale  = {float(robin_ratio(36, mpf(2))):.10f}")

        section("2. Numerical critical exponent of n = 36")
        ks = kappa_star(36)
        print(f"  kappa*(36) = {float(ks):.6f}")
        print(f"  ratio below the root = {float(robin_ratio(36, ks - mpf('0.01'))):.6f}")
        print(f"  ratio above the root = {float(robin_ratio(36, ks + mpf('0.01'))):.6f}")

        section("3. Finite exceptional-set search at kappa = 2")
        exceptions = exceptional_set(mpf(2), n_max=100)
        print(f"  E(2) intersect [3, 100] = {exceptions}")

        section("4. A reciprocal 2-CA number at epsilon = -0.40")
        fac: Factorization = reciprocal_ca(mpf(2), mpf("-0.40"))
        print(f"  omega                      = {fac.omega}")
        print(f"  largest prime divisor      = {fac.p_max}")
        print(f"  log n                      = {float(fac.log_n()):.6f}")
        print(f"  rho^[2] / Robin scale at N = {float(fac.robin_ratio(mpf(2))):.6f}")
        print("  This is a pointwise CA computation; it does not control nearby integers.")

    print("\nDone. See ARCHITECTURE.md and CERTIFICATION.md for scope.")


if __name__ == "__main__":
    main()
