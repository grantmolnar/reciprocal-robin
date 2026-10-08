"""Generate and verify certificates at progressively larger integer ``K``.

For each parameter, the sweep generates a fresh certificate, runs the independent
verifier, and only then records certificate statistics and advances. A timeout or
failed process marks the limit of that computation; it is not a mathematical
counterexample.

A source checkout invokes the wrappers under ``scripts/``. An installed package
invokes the corresponding modules with ``python -m``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shlex
import subprocess
import sys
import time
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


@dataclass(frozen=True)
class StageResult:
    status: str
    seconds: float
    returncode: int | None
    log_path: Path


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_values(text: str) -> list[int]:
    values: list[int] = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            value = int(item)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                f"every K must be an integer > 1; got {item!r}"
            ) from exc
        if value <= 1:
            raise argparse.ArgumentTypeError("every K must be an integer > 1")
        values.append(value)
    if not values:
        raise argparse.ArgumentTypeError("--values must contain at least one integer")
    if any(b <= a for a, b in zip(values, values[1:])):
        raise argparse.ArgumentTypeError("--values must be strictly increasing")
    return values


def geometric_values(start: int, factor: float, max_k: int | None) -> Iterable[int]:
    k = start
    while max_k is None or k <= max_k:
        yield k
        next_k = int(math.ceil(k * factor))
        if next_k <= k:
            next_k = k + 1
        k = next_k


def run_stage(
    command: list[str],
    *,
    cwd: Path,
    log_path: Path,
    timeout_seconds: float | None,
) -> StageResult:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"started_utc: {utc_now()}\n")
        log.write("command: " + shlex.join(command) + "\n\n")
        log.flush()
        try:
            completed = subprocess.run(
                command,
                cwd=cwd,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=timeout_seconds,
                check=False,
            )
        except subprocess.TimeoutExpired:
            seconds = time.perf_counter() - started
            log.write(f"\nTIMEOUT after {seconds:.3f} seconds\n")
            return StageResult("timeout", seconds, None, log_path)
        except OSError as exc:
            seconds = time.perf_counter() - started
            log.write(f"\nOS ERROR after {seconds:.3f} seconds: {exc!r}\n")
            return StageResult("os_error", seconds, None, log_path)

    seconds = time.perf_counter() - started
    status = "ok" if completed.returncode == 0 else "failed"
    return StageResult(status, seconds, completed.returncode, log_path)


def read_certificate_stats(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    tail = data["finite_tail"]
    endpoints = data["endpoint_verification"]
    initial = data["initial_integer_scan"]
    return {
        "tail_start_prime": tail["tail_start_prime"],
        "last_bad_prime": tail["last_bad_prime"],
        "prime_rows_checked": tail["prime_rows_checked"],
        "canonical_events": endpoints["event_count"],
        "first_endpoint_above_5040": endpoints["first_endpoint_above_5040"],
        "minimum_endpoint_log_margin_lo": endpoints["minimum_log_margin_lo"],
        "minimum_endpoint_log_margin_endpoint": endpoints["minimum_log_margin_endpoint"],
        "initial_integer_checks": initial["count"],
        "initial_minimum_margin_lo": initial["minimum_margin_lo"],
        "initial_minimum_at_n": initial["minimum_at_n"],
    }


def append_csv(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "timestamp_utc",
        "kappa",
        "status",
        "dps",
        "timeout_minutes",
        "generate_seconds",
        "verify_seconds",
        "certificate_bytes",
        "certificate_sha256",
        "tail_start_prime",
        "last_bad_prime",
        "prime_rows_checked",
        "canonical_events",
        "first_endpoint_above_5040",
        "minimum_endpoint_log_margin_lo",
        "minimum_endpoint_log_margin_endpoint",
        "initial_integer_checks",
        "initial_minimum_margin_lo",
        "initial_minimum_at_n",
        "generator_log",
        "verifier_log",
    ]
    exists = path.exists()
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def append_jsonl(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        json.dump(row, handle, sort_keys=True)
        handle.write("\n")


def _proof_commands(python: str, repo: Path | None) -> tuple[list[str], list[str], Path]:
    if repo is not None:
        generator = repo / "scripts" / "certify_bounded_parameter_interval.py"
        verifier = repo / "scripts" / "verify_bounded_parameter_certificate.py"
        if not generator.is_file() or not verifier.is_file():
            raise ValueError("repository root does not contain the bounded-parameter proof scripts")
        return [python, str(generator)], [python, str(verifier)], repo
    return (
        [python, "-m", "reciprocal_robin.certification.bounded_generate"],
        [python, "-m", "reciprocal_robin.certification.bounded_verify"],
        Path.cwd(),
    )


def _verified_values(csv_path: Path) -> set[int]:
    if not csv_path.is_file():
        return set()
    values: set[int] = set()
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("status") == "verified":
                values.add(int(row["kappa"]))
    return values


def main(argv: Sequence[str] | None = None, *, default_repo: Path | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Generate and independently verify bounded-parameter certificates "
            "at progressively larger integer K."
        )
    )
    parser.add_argument(
        "--repo",
        type=Path,
        default=default_repo,
        help=(
            "source-checkout repository root; omit when using an installed package "
            "(the source wrapper supplies this automatically)"
        ),
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--values",
        type=parse_values,
        help="explicit increasing comma-separated K values, e.g. 2000,5000,10000",
    )
    group.add_argument(
        "--start",
        type=int,
        default=2000,
        help="first K for geometric sweep (default: 2000)",
    )
    parser.add_argument(
        "--factor",
        type=float,
        default=2.0,
        help="geometric growth factor when --values is not used (default: 2.0)",
    )
    parser.add_argument(
        "--max-k",
        type=int,
        default=1_000_000,
        help="largest K attempted in geometric mode (default: 1000000)",
    )
    parser.add_argument(
        "--dps",
        type=int,
        default=60,
        help="directed-interval decimal precision passed to generator (default: 60)",
    )
    parser.add_argument(
        "--timeout-minutes",
        type=float,
        default=30.0,
        help="timeout for EACH generator/verifier stage; 0 means no timeout (default: 30)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("sweep-results"),
        help="directory for certificates, logs, and ledgers (default: sweep-results)",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="skip K values already recorded as verified in sweep.csv",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python interpreter used for proof scripts (default: current interpreter)",
    )
    args = parser.parse_args(argv)

    repo = args.repo.resolve() if args.repo is not None else None
    try:
        generator_command, verifier_command, cwd = _proof_commands(args.python, repo)
    except ValueError as exc:
        parser.error(str(exc))
    if args.dps < 2:
        parser.error("--dps must be at least 2")
    if not math.isfinite(args.timeout_minutes) or args.timeout_minutes < 0:
        parser.error("--timeout-minutes must be finite and nonnegative")

    if args.values is None:
        if args.start <= 1:
            parser.error("--start must be > 1")
        if not math.isfinite(args.factor) or args.factor <= 1.0:
            parser.error("--factor must be finite and > 1")
        if args.max_k is not None and args.max_k < args.start:
            parser.error("--max-k must be >= --start")
        values: Iterable[int] = geometric_values(args.start, args.factor, args.max_k)
    else:
        values = args.values

    output_dir = args.output_dir
    if not output_dir.is_absolute():
        output_dir = cwd / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = output_dir / "logs"
    certs_dir = output_dir / "certificates"
    csv_path = output_dir / "sweep.csv"
    jsonl_path = output_dir / "sweep.jsonl"
    timeout_seconds = None if args.timeout_minutes == 0 else 60.0 * args.timeout_minutes

    already_verified = _verified_values(csv_path) if args.resume else set()

    execution = f"source checkout {repo}" if repo is not None else "installed package"
    print(f"Execution:  {execution}")
    print(f"Python:     {args.python}")
    print(f"Results:    {output_dir}")
    print(f"DPS:        {args.dps}")
    print(
        "Timeout:    "
        + ("none" if timeout_seconds is None else f"{args.timeout_minutes:g} min per stage")
    )
    print()

    last_success: int | None = max(already_verified, default=None)
    for kappa in values:
        if kappa in already_verified:
            print(f"=== K={kappa}: already verified; skipping (--resume) ===")
            last_success = max(last_success or kappa, kappa)
            continue
        certificate = certs_dir / f"kappa_{kappa}_endpoint_certificate.json"
        generator_log = logs_dir / f"kappa_{kappa}_generate.log"
        verifier_log = logs_dir / f"kappa_{kappa}_verify.log"

        print(f"=== K={kappa}: generating certificate ===", flush=True)
        generate = run_stage(
            [
                *generator_command,
                "--kappa",
                str(kappa),
                "--dps",
                str(args.dps),
                "--output",
                str(certificate),
            ],
            cwd=cwd,
            log_path=generator_log,
            timeout_seconds=timeout_seconds,
        )

        base_row: dict[str, object] = {
            "timestamp_utc": utc_now(),
            "kappa": kappa,
            "dps": args.dps,
            "timeout_minutes": args.timeout_minutes,
            "generate_seconds": round(generate.seconds, 6),
            "verify_seconds": "",
            "certificate_bytes": "",
            "certificate_sha256": "",
            "generator_log": os.path.relpath(generator_log, cwd),
            "verifier_log": os.path.relpath(verifier_log, cwd),
        }

        if generate.status != "ok":
            row = {**base_row, "status": f"generator_{generate.status}"}
            append_csv(csv_path, row)
            append_jsonl(jsonl_path, row)
            print(
                f"STOP at K={kappa}: generator {generate.status} " f"after {generate.seconds:.1f}s."
            )
            print(f"Log: {generator_log}")
            if last_success is not None:
                print(f"Largest independently verified K in this sweep: {last_success}")
            return 2

        try:
            cert_size = certificate.stat().st_size
            cert_sha = sha256_file(certificate)
        except OSError as exc:
            row = {**base_row, "status": "generator_missing_output"}
            append_csv(csv_path, row)
            append_jsonl(jsonl_path, row)
            print(f"STOP at K={kappa}: cannot read generated certificate: {exc}")
            print(f"Log: {generator_log}")
            return 2

        print(
            f"Generated in {generate.seconds:.1f}s; "
            f"{cert_size / (1024 * 1024):.2f} MiB; verifying...",
            flush=True,
        )

        verify = run_stage(
            [
                *verifier_command,
                "--certificate",
                str(certificate),
            ],
            cwd=cwd,
            log_path=verifier_log,
            timeout_seconds=timeout_seconds,
        )

        row = {
            **base_row,
            "verify_seconds": round(verify.seconds, 6),
            "certificate_bytes": cert_size,
            "certificate_sha256": cert_sha,
        }
        if verify.status != "ok":
            row["status"] = f"verifier_{verify.status}"
            append_csv(csv_path, row)
            append_jsonl(jsonl_path, row)
            print(f"STOP at K={kappa}: verifier {verify.status} " f"after {verify.seconds:.1f}s.")
            print(f"Log: {verifier_log}")
            if last_success is not None:
                print(f"Largest independently verified K in this sweep: {last_success}")
            return 3

        # Read telemetry only after verification. Advancement depends on the
        # verifier, not on this summary code's reading of the certificate.
        try:
            stats = read_certificate_stats(certificate)
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
            stats = {}
            print(f"Warning: verified certificate statistics could not be summarized: {exc}")

        row.update(stats)
        row["status"] = "verified"
        append_csv(csv_path, row)
        append_jsonl(jsonl_path, row)
        last_success = kappa
        summary = (
            f"tail P={stats['tail_start_prime']}, events={stats['canonical_events']}, "
            f"initial checks={stats['initial_integer_checks']}"
            if stats
            else "certificate statistics unavailable"
        )
        print(
            f"VERIFIED K={kappa} in {generate.seconds + verify.seconds:.1f}s total; " f"{summary}."
        )
        print()

    print("Sweep completed without hitting the configured frontier.")
    if last_success is not None:
        print(f"Largest independently verified K: {last_success}")
    print(f"Ledger: {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
