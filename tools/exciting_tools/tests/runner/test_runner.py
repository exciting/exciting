"""Tests for the binary runner.
Excludes the run method.
"""
import shutil
from pathlib import Path
from typing import Any

import pytest

from excitingtools.runner.runner import BinaryRunner, RunnerCode

mock_binary = "false_exciting_binary"


@pytest.mark.xfail(shutil.which(mock_binary) is not None, reason="Binary name exists.")
def test_no_binary():
    with pytest.raises(FileNotFoundError,
                       match=fr"{mock_binary} binary is not present in the current directory nor in \$PATH"):
        BinaryRunner(mock_binary, "./", 1, 1)


@pytest.fixture
def exciting_smp(tmp_path: Path) -> str:
    binary = tmp_path / "exciting_smp"
    binary.touch()
    return binary.as_posix()


def test_no_run_dir(exciting_smp: str):
    with pytest.raises(OSError, match="Run directory does not exist: non_existent_dir"):
        BinaryRunner(exciting_smp, "./", 1, 1, "non_existent_dir")


def test_false_run_cmd(exciting_smp: str):
    false_run_cmd: Any = 3
    with pytest.raises(ValueError, match="Run commands expected in a str or list. For example ['mpirun', '-np', '2']"):
        BinaryRunner(exciting_smp, false_run_cmd, 1, 1)


@pytest.fixture
def exciting_mpismp(tmp_path: Path) -> str:
    binary = tmp_path / "exciting_mpismp"
    binary.touch()
    return binary.as_posix()


def test_false_mpi_command_smaller_than_zero(exciting_mpismp: str):
    with pytest.raises(ValueError, match="Number of MPI processes must be > 0"):
        BinaryRunner(exciting_mpismp, ["mpirun", "-np", "-1"], 1, 1)


def test_false_mpi_command_no_int(exciting_mpismp: str):
    with pytest.raises(ValueError, match="Number of MPI processes should be an int"):
        BinaryRunner(exciting_mpismp, ["mpirun", "-np", "no_int"], 1, 1)


def test_false_mpi_command_no_number_given(exciting_mpismp: str):
    with pytest.raises(ValueError, match="Number of MPI processes must be specified after the '-np'"):
        BinaryRunner(exciting_mpismp, ["mpirun", "-np"], 1, 1)


def test_false_omp(exciting_smp: str):
    with pytest.raises(ValueError, match="Number of OMP threads must be > 0"):
        BinaryRunner(exciting_smp, [""], -1, 1)


def test_false_timeout(exciting_smp: str):
    with pytest.raises(ValueError, match="time_out must be a positive integer"):
        BinaryRunner(exciting_smp, [""], 1, -1)


@pytest.fixture
def runner(tmp_path: Path, exciting_mpismp: str) -> BinaryRunner:
    """Produces a runner with binary and run dir mocked up.
    """
    run_dir = tmp_path / "ab/de"
    run_dir.mkdir(parents=True)
    return BinaryRunner(exciting_mpismp, ["mpirun", "-np", "3"], 4, 260, run_dir.as_posix(), [">", "std.out"])


def test_as_dict(tmp_path: Path, runner: BinaryRunner, mock_env_jobflow_missing):
    assert runner.as_dict() == {'args': ['>', 'std.out'],
                                'binary': (tmp_path / "exciting_mpismp").as_posix(),
                                'directory': (tmp_path / "ab/de").as_posix(),
                                'omp_num_threads': 4,
                                'run_cmd': ['mpirun', '-np', '3'],
                                'time_out': 260}


def test_as_dict_jobflow(tmp_path: Path, runner: BinaryRunner, mock_env_jobflow):
    assert runner.as_dict() == {'@class': 'BinaryRunner',
                                '@module': 'excitingtools.runner.runner',
                                'args': ['>', 'std.out'],
                                'binary': (tmp_path / "exciting_mpismp").as_posix(),
                                'directory': (tmp_path / "ab/de").as_posix(),
                                'omp_num_threads': 4,
                                'run_cmd': ['mpirun', '-np', '3'],
                                'time_out': 260}


def test_from_dict(tmp_path: Path, runner):
    new_runner = BinaryRunner.from_dict(runner.as_dict())
    assert new_runner.binary == (tmp_path / "exciting_mpismp").as_posix()
    assert new_runner.time_out == 260


def test_run_with_bash_command(tmp_path: Path):
    """Produces a runner with binary and run dir mocked up.
    Test a simple echo command.
    """
    run_dir = tmp_path / "ab/de"
    run_dir.mkdir(parents=True)
    binary = tmp_path / "hello"
    binary.touch()
    runner = BinaryRunner(binary.as_posix(), ["echo"], 1, 60, run_dir.as_posix())
    run_results = runner.run()
    assert run_results.success
    assert run_results.stderr == ''
    assert run_results.stdout == binary.as_posix() + "\n"


def test_timeout_with_bash_command(tmp_path: Path):
    """Produces a runner with binary and run dir mocked up.

    Test a simple sleep command to get a timeout.
    """
    time_out = 1
    binary = tmp_path / "sleep.sh"
    binary.write_text(f"sleep {time_out + 0.1}")
    runner = BinaryRunner(binary.as_posix(), ["sh"], 1, time_out, tmp_path.as_posix())
    run_results = runner.run()
    assert not run_results.success
    assert run_results.stderr == 'BinaryRunner: Job timed out. \n\n'
    assert run_results.stdout == ""
    assert run_results.process_time == time_out
    assert isinstance(run_results.return_code, RunnerCode)
    assert run_results.return_code == RunnerCode.time_out
