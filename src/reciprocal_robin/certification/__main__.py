"""Command line for certificate generation and verification.

Dispatch imports only the selected command module. A verifier invocation
therefore does not import its corresponding generator through this entry point.
"""

from __future__ import annotations

import sys
from collections.abc import Callable, Sequence
from importlib import import_module
from typing import cast

_Command = Callable[[Sequence[str] | None], int]

_COMMANDS: dict[str, tuple[str, str]] = {
    "generate": ("bounded_generate", "generate one directed-interval certificate"),
    "verify": ("bounded_verify", "independently verify one certificate"),
    "sweep": ("sweep", "generate and verify progressively larger integer K"),
    "generate-kappa-five": (
        "kappa_five_generate",
        "generate the shared kappa<=5 finite certificate",
    ),
    "verify-kappa-five": (
        "kappa_five_verify",
        "independently verify that certificate",
    ),
}


def _usage() -> str:
    names = ",".join(_COMMANDS)
    rows = "\n".join(f"  {name:<20} {description}" for name, (_, description) in _COMMANDS.items())
    return (
        "usage: python -m reciprocal_robin.certification "
        f"{{{names}}} [options]\n\ncommands:\n{rows}\n\n"
        "Run a command with --help for its options."
    )


def _load_command(module_name: str) -> _Command:
    module = import_module(f"reciprocal_robin.certification.{module_name}")
    return cast(_Command, getattr(module, "main"))


def main(argv: Sequence[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if not arguments or arguments[0] in {"-h", "--help"}:
        print(_usage())
        return 0
    command = arguments.pop(0)
    command_spec = _COMMANDS.get(command)
    if command_spec is None:
        print(f"unknown certification command: {command!r}\n\n{_usage()}", file=sys.stderr)
        return 2
    module_name, _ = command_spec
    return _load_command(module_name)(arguments)


if __name__ == "__main__":
    raise SystemExit(main())
