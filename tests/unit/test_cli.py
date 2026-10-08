"""In-process tests for the CLI subcommands.

We call ``cli.main([...])`` directly rather than spawning a subprocess.
This is faster, gives coverage instrumentation a chance to track the
calls, and avoids depending on the install location of the entry point.
"""

import pytest

from reciprocal_robin.cli import main


def test_cli_pointwise(capsys):
    rc = main(["pointwise", "36", "2"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "0.992821" in out


def test_cli_critical_exponent(capsys):
    rc = main(["critical-exponent", "36"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "7.81" in out


def test_cli_exceptional_set(capsys):
    rc = main(["exceptional-set", "2", "--max-n", "30"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "[3, 4, 5, 6, 8, 12, 16, 24]" in out


def test_cli_walk(capsys):
    rc = main(
        [
            "walk",
            "2",
            "--eps-min",
            "-0.40",
            "--eps-max",
            "-0.30",
            "--steps",
            "3",
        ]
    )
    assert rc == 0
    out = capsys.readouterr().out
    assert "omega=" in out


def test_cli_scan_grid_returns_zero_when_sampled_ratios_are_below_one(capsys):
    rc = main(
        [
            "scan-grid",
            "2",
            "--eps-min",
            "-0.40",
            "--eps-max",
            "-0.30",
            "--steps",
            "3",
        ]
    )
    assert rc == 0
    out = capsys.readouterr().out
    assert "Max ratio observed" in out


def test_cli_unknown_subcommand_exits_nonzero():
    with pytest.raises(SystemExit) as exc_info:
        main(["nonexistent-subcommand"])
    assert exc_info.value.code != 0


@pytest.mark.parametrize(
    "argv, message",
    [
        (["pointwise", "2", "2"], "n must be at least 3"),
        (["pointwise", "36", "1"], "kappa must be greater than"),
        (["walk", "2", "--eps-min", "nan", "--eps-max", "-0.3"], "eps-min must be finite"),
        (
            ["scan-grid", "2", "--eps-min", "-0.4", "--eps-max", "-0.3", "--steps", "0"],
            "steps must be at least 1",
        ),
    ],
)
def test_cli_rejects_invalid_numeric_arguments_cleanly(argv, message, capsys):
    with pytest.raises(SystemExit) as exc_info:
        main(argv)
    assert exc_info.value.code == 2
    err = capsys.readouterr().err
    assert message in err
    assert "Traceback" not in err


def test_cli_rejects_epsilon_outside_ca_domain_cleanly(capsys):
    with pytest.raises(SystemExit) as exc_info:
        main(
            [
                "walk",
                "2",
                "--eps-min",
                "-0.6",
                "--eps-max",
                "-0.55",
                "--steps",
                "2",
            ]
        )
    assert exc_info.value.code == 2
    err = capsys.readouterr().err
    assert "epsilon must be greater than -1/kappa" in err
    assert "Traceback" not in err
