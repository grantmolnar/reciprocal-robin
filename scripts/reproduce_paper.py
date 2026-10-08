"""Check every registered numerical claim from the manuscript.

The script runs ``tests/paper/test_paper_numerics.py`` and returns nonzero if any
registered value no longer agrees with the code.
"""

import argparse
import os
import subprocess
import sys
from collections.abc import Sequence

from release_support import REPO_ROOT, SRC_DIR

TEST_FILES = [
    REPO_ROOT / "tests" / "paper" / "test_paper_numerics.py",
    REPO_ROOT / "tests" / "certification" / "test_kappa_five_interval.py",
    REPO_ROOT / "tests" / "certification" / "test_certification_results.py",
]


def _build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        description=(
            "Reproduce the registered manuscript numerics and compact "
            "certificate checks. Unrecognized options are passed "
            "through to pytest."
        )
    )


def main(argv: Sequence[str] | None = None) -> int:
    _, pytest_args = _build_parser().parse_known_args(argv)
    missing = [path for path in TEST_FILES if not path.exists()]
    if missing:
        print(f"FATAL: missing test files: {missing}", file=sys.stderr)
        return 2
    extra = pytest_args or ["-v"]
    print("Reproducing manuscript numerics and compact certificates ...")
    env = os.environ.copy()
    src = str(SRC_DIR)
    env["PYTHONPATH"] = src if not env.get("PYTHONPATH") else src + os.pathsep + env["PYTHONPATH"]
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", *(str(path) for path in TEST_FILES), *extra],
        cwd=REPO_ROOT,
        env=env,
    )
    if proc.returncode == 0:
        print("\nAll registered manuscript numerics reproduce.")
    else:
        print(
            "\nDRIFT DETECTED: at least one numerical claim in the paper no "
            "longer matches the code. Either the code regressed, or a paper "
            "revision changed a printed constant without updating the test.",
            file=sys.stderr,
        )
    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
