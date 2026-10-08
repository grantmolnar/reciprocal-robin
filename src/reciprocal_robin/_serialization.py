"""Stable UTF-8 JSON serialization shared by reports and certificates.

This module contains no theorem-specific arithmetic. Independent verifiers may
therefore import it without depending on certificate generators.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any


def json_text(payload: Any, *, trailing_newline: bool = False) -> str:
    """Serialize ``payload`` in the repository's stable human-readable format."""
    text = json.dumps(payload, indent=2, ensure_ascii=False)
    return text + "\n" if trailing_newline else text


def read_json(path: Path) -> Any:
    """Read a UTF-8 JSON document from ``path``."""
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    """Atomically write a UTF-8 JSON document with one trailing newline."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            handle.write(json_text(payload, trailing_newline=True))
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
