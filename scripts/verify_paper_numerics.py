"""Evaluate every numerical claim registered for the manuscript.

Each claim receives a one-line verdict. The command succeeds only when all
registered claims are confirmed.
"""

import argparse
import sys
from collections.abc import Sequence

from release_support import REPO_ROOT  # noqa: F401
from reciprocal_robin.paper_numerics import Status, verify_all


def _build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        description="Verify every numerical claim registered from the companion paper."
    )


def main(argv: Sequence[str] | None = None) -> int:
    _build_parser().parse_args(argv)
    verdicts = verify_all()
    width = max(len(v.claim.locator) for v in verdicts)
    counts = {s: 0 for s in Status}

    for v in verdicts:
        counts[v.status] += 1
        tag = {
            Status.CONFIRMED: "OK   ",
            Status.SLIP: "SLIP ",
            Status.UNKNOWN: "???? ",
        }[v.status]
        print(f"[{tag}] {v.claim.locator:<{width}}  {v.claim.description}")
        if v.note:
            print(f"         {v.note}")

    print()
    print(
        f"Summary: {counts[Status.CONFIRMED]} confirmed, "
        f"{counts[Status.SLIP]} slip, {counts[Status.UNKNOWN]} unknown "
        f"({len(verdicts)} total)"
    )
    return 0 if counts[Status.SLIP] == 0 and counts[Status.UNKNOWN] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
