"""Tests for the monitored runner."""

import re
import sys
from pathlib import Path

import pytest

from excitingtools.runner.monitored_runner import MonitoredGroundStateRunner, RunnerCode


@pytest.fixture
def python_file(tmp_path: Path) -> Path:
    file = tmp_path / "test.py"
    file.write_text(r"""from time import sleep
with open("INFO.OUT", "w+") as f:
    sleep(2)
    f.write("SCF iteration number : 1\n\n\n")
    f.write("SCF iteration number : 2\n\n\n")
    f.write("SCF iteration number : 3\n\n\n")
    f.write("stopped")
""")
    return file


def test_runner(tmp_path: Path, capsys, python_file: Path):
    """Tests the monitored runner and its output."""
    runner = MonitoredGroundStateRunner(python_file, sys.executable, directory=tmp_path)

    run_results = runner.run()

    # get all the printed lines from the captured sysout
    lines = capsys.readouterr().out.splitlines()
    """The expected lines looks like (numbers are examples and may vary between runs):
        lines = ["Exciting started!",
        "Iteration 1 (0.0s)",
        "\033[FIteration 1 (0.1s)",
        "\033[FIteration 1 (0.2s)",
        ...
        "\033[FLast finished iteration:   1 (2.9s)",
        "Iteration 2 (0.0s)",
        "\033[F\033[FLast finished iteration:   2 (0.1s)",
        "Iteration 3 (0.0s)",
        "\033[F\033[FLast finished iteration:   3 (0.1s)",
        "Iteration 4 (0.0s)",
        "Exciting Finished! (0m3s)"]
        
        Note: we only have to wait for the first iteration, because we can not access the file, but as soon as the file
        is readable we can read everything and the program will rapidly finish reading.
    """

    # regex pattern for a line which gets printed for a currently running iteration and for a finished iteration
    running_iteration = re.compile(r"(?:\033\[F)?Iteration (\d+) \(\d+\.\ds\)")
    finished_iteration = re.compile(r"(?:\033\[F){1,2}Last finished iteration:\s+(\d+) \(\d+\.\ds\)")

    current_iteration = 1
    assert lines[0] == "Exciting started!"
    # iterate over all lines, which are printed while the program is running
    for line in lines[1:-1]:
        match = running_iteration.match(line) or finished_iteration.match(line)
        assert match is not None, f"{line!r} does not match expected output!"
        assert current_iteration == int(match.group(1)), (
            "Iteration number mismatch! " f"Expected {current_iteration}, got {int(match.group(1))}"
        )
        if "finished" in line:
            current_iteration += 1

    assert re.match(r"Exciting Finished! \(\d+m\d\ds\)", lines[-1]) is not None

    assert run_results.success


def test_timeout(tmp_path: Path, capsys, python_file: Path):
    """Test the monitored runner to get a timeout."""
    runner = MonitoredGroundStateRunner(python_file, sys.executable, directory=tmp_path, time_out=1)

    run_results = runner.run()

    lines = capsys.readouterr().out.splitlines()
    assert lines[-2] == "Exciting timed out!"

    assert isinstance(run_results.return_code, RunnerCode)
    assert run_results.return_code == RunnerCode.time_out
