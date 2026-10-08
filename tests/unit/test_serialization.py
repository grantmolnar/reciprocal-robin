"""Serialization contracts shared by certificates and release metadata."""

from pathlib import Path

from reciprocal_robin._serialization import read_json, write_json


def test_json_helpers_round_trip_with_trailing_newline(tmp_path: Path) -> None:
    path = tmp_path / "certificate.json"
    payload = {"schema_version": "test-1.0", "rows": [{"n": 3}]}
    write_json(path, payload)
    assert path.read_bytes().endswith(b"\n")
    assert read_json(path) == payload
    assert not list(tmp_path.glob("*.tmp"))
