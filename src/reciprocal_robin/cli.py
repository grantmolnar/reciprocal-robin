"""Command-line interface for pointwise and sampled calculations.

Certificate generation and verification use the separate entry points described
in :mod:`reciprocal_robin.certification` and ``CERTIFICATION.md``.

Examples
--------
::

    reciprocal-robin pointwise 36 2
    reciprocal-robin critical-exponent 48
    reciprocal-robin exceptional-set 2 --max-n 100
    reciprocal-robin walk 2 --eps-min -0.49 --eps-max -0.30 --steps 10
    reciprocal-robin scan-grid 2 --eps-min -0.49 --eps-max -0.30 --steps 10
    reciprocal-robin audit-grid 2 --eps-min -0.49 --eps-max -0.30 --steps 10
"""

import argparse
import sys
from collections.abc import Callable, Sequence
from dataclasses import asdict
from pathlib import Path

import mpmath
from mpmath import mpf

from . import MpfField, exceptional_set, kappa_star, rho_k, robin_ratio, using_field
from ._precision import DEFAULT_DPS, working_precision
from ._serialization import write_json
from ._types import Scalar
from ._validation import require_integer, require_kappa
from .exploration import (
    audit_ca_grid,
    max_robin_ratio_on_ca_grid,
    print_step,
    reciprocal_ca,
)


def _integer_at_least(name: str, minimum: int) -> Callable[[str], int]:
    """Return an argparse converter for an integer with a lower bound."""

    def parse(text: str) -> int:
        try:
            value = int(text)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"{name} must be an integer") from exc
        if value < minimum:
            raise argparse.ArgumentTypeError(f"{name} must be at least {minimum}")
        return value

    return parse


def _kappa_text(text: str) -> str:
    """Validate a CLI kappa while retaining the user's text for display."""
    try:
        value = mpf(text)
        require_kappa(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return text


def _finite_real_text(name: str) -> Callable[[str], str]:
    """Return an argparse converter for a finite real while preserving text."""

    def parse(text: str) -> str:
        try:
            value = mpf(text)
        except (TypeError, ValueError) as exc:
            raise argparse.ArgumentTypeError(f"{name} must be a real number") from exc
        if not mpmath.isfinite(value):
            raise argparse.ArgumentTypeError(f"{name} must be finite")
        return text

    return parse


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="reciprocal-robin",
        description="Numerical exploration of the reciprocal sigma family.",
    )
    subparsers = parser.add_subparsers(dest="cmd", required=True)

    command = subparsers.add_parser("pointwise", help="rho^[k](n) and the Robin ratio")
    command.add_argument("n", type=_integer_at_least("n", 3))
    command.add_argument("k", type=_kappa_text, help="exponent kappa > 1, e.g. 2 or 1.5")

    command = subparsers.add_parser("critical-exponent", help="kappa*(n)")
    command.add_argument("n", type=_integer_at_least("n", 3))

    command = subparsers.add_parser(
        "exceptional-set",
        help="integers in [3, n_max] violating the inequality",
    )
    command.add_argument("k", type=_kappa_text)
    command.add_argument("--max-n", type=_integer_at_least("max-n", 3), required=True)

    command = subparsers.add_parser(
        "walk",
        help="enumerate reciprocal kappa-CA numbers across an epsilon grid",
    )
    command.add_argument("k", type=_kappa_text)
    command.add_argument("--eps-min", type=_finite_real_text("eps-min"), required=True)
    command.add_argument("--eps-max", type=_finite_real_text("eps-max"), required=True)
    command.add_argument("--steps", type=_integer_at_least("steps", 1), default=10)

    command = subparsers.add_parser(
        "scan-grid",
        help="report the maximum point ratio on an explicit grid",
    )
    command.add_argument("k", type=_kappa_text)
    command.add_argument("--eps-min", type=_finite_real_text("eps-min"), required=True)
    command.add_argument("--eps-max", type=_finite_real_text("eps-max"), required=True)
    command.add_argument("--steps", type=_integer_at_least("steps", 1), default=10)

    command = subparsers.add_parser(
        "audit-grid",
        help="emit interval verdicts for sampled reciprocal-CA points (JSON)",
    )
    command.add_argument("k", type=_kappa_text)
    command.add_argument("--eps-min", type=_finite_real_text("eps-min"), required=True)
    command.add_argument("--eps-max", type=_finite_real_text("eps-max"), required=True)
    command.add_argument("--steps", type=_integer_at_least("steps", 1), default=10)
    command.add_argument(
        "--base-dps",
        type=_integer_at_least("base-dps", 2),
        default=50,
    )
    command.add_argument(
        "--bumped-dps",
        type=_integer_at_least("bumped-dps", 2),
        default=200,
    )
    command.add_argument(
        "--output",
        type=Path,
        default=None,
        help="write JSON to this file (default: stdout)",
    )

    return parser


def _eps_grid(lo: str, hi: str, steps: int) -> list[Scalar]:
    """Return an inclusive, linearly spaced epsilon grid."""
    require_integer(steps, "steps", minimum=1)
    lower, upper = mpf(lo), mpf(hi)
    if steps == 1:
        return [lower]
    return [lower + (upper - lower) * mpf(index) / (steps - 1) for index in range(steps)]


def _main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        if args.cmd == "pointwise":
            kappa = mpf(args.k)
            with using_field(MpfField()):
                rho = rho_k(args.n, kappa)
                ratio = robin_ratio(args.n, kappa)
            print(f"rho^[{args.k}]({args.n})        = {float(rho):.12f}")
            print(f"robin_ratio({args.n}, {args.k}) = {float(ratio):.12f}")
            return 0

        if args.cmd == "critical-exponent":
            with using_field(MpfField()):
                critical = kappa_star(args.n)
            print(f"kappa*({args.n}) = {float(critical):.6f}")
            return 0

        if args.cmd == "exceptional-set":
            with using_field(MpfField()):
                exceptions = exceptional_set(mpf(args.k), n_max=args.max_n)
            print(f"E({args.k}) cap [3, {args.max_n}] = {exceptions}")
            return 0

        if args.cmd == "walk":
            grid = _eps_grid(args.eps_min, args.eps_max, args.steps)
            with using_field(MpfField()):
                rows = []
                for eps in grid:
                    factorization = reciprocal_ca(mpf(args.k), eps)
                    rows.append(
                        (
                            eps,
                            factorization.omega,
                            factorization.p_max,
                            factorization.log_n(),
                        )
                    )
            for eps, omega, p_max, log_n in rows:
                print(
                    f"eps={float(eps):+.4f}  omega={omega:6d}  P*={p_max:10d}  "
                    f"log n={float(log_n):12.2f}"
                )
            return 0

        if args.cmd == "scan-grid":
            grid = _eps_grid(args.eps_min, args.eps_max, args.steps)
            with using_field(MpfField()):
                max_ratio = max_robin_ratio_on_ca_grid(args.k, grid, observer=print_step)
            print(f"\nMax ratio observed: {float(max_ratio):.6f}")
            # This exit code concerns only the sampled points.
            return 0 if float(max_ratio) < 1 else 1

        if args.cmd == "audit-grid":
            # Parse the grid at the same precision used for the interval audit;
            # otherwise high-precision evaluation cannot recover endpoint digits.
            with working_precision(args.base_dps):
                grid = _eps_grid(args.eps_min, args.eps_max, args.steps)
            report = audit_ca_grid(
                args.k,
                grid,
                base_dps=args.base_dps,
                bumped_dps=args.bumped_dps,
            )
            payload = report.to_json()
            if args.output is not None:
                write_json(args.output, asdict(report))
                print(
                    f"wrote {args.output}: kappa={report.kappa}, "
                    f"steps={len(report.steps)}, "
                    f"resolved={report.all_steps_resolved}",
                    file=sys.stderr,
                )
            else:
                print(payload)
            return 0 if report.all_steps_resolved else 1
    except (TypeError, ValueError) as exc:
        parser.error(str(exc))

    return 2  # argparse requires one of the registered subcommands


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI at the package baseline precision."""
    with working_precision(DEFAULT_DPS):
        return _main(argv)


if __name__ == "__main__":
    sys.exit(main())
