"""Monitored runner, printer for multiple changing lines and asynchronous functions for real time monitoring of the
INFO.OUT file."""

import asyncio
import re
import subprocess
from datetime import timedelta
from pathlib import Path
from time import perf_counter
from typing import List

from excitingtools.runner.runner import BinaryRunner, RunnerCode, SubprocessRunResults


class ExcitingException(BaseException):
    """
    Exception class if the exciting process stops unexpectedly.
    """


class Sentinel:
    """
    Class used for signaling the end of the async I/O process.
    """


class MultiLinePrinter:
    """A class used to print multiple lines of text where each subsequent print overwrites the previous batch and
    completely overwrites the previous text even if the new text is smaller than the previous one.
    """

    def __init__(self):
        """
        Initializes a MultiLinePrinter object.
        """
        self.line_limits = []

    def print(self, lines: str):
        """
        Prints the input lines in a formatted manner. If a line is shorter than its previous print, it will be padded
        with spaces. If there are fewer lines than before, additional blank lines will be printed. The padding ensures
        that no leftover text from the previous print remains.

        :param lines: The lines of text to be printed. The lines should be separated by newline characters.
        """
        # Move the cursor up for the number of lines previously printed
        print("\033[F" * len(self.line_limits), end="")  # noqa: T201

        # Split the input string into separate lines
        lines = lines.splitlines()
        for i, line in enumerate(lines):
            # If this is a new line, add its length to the line limits
            if i == len(self.line_limits):
                self.line_limits.append(len(line))
            else:
                # If this line is longer than before, update the line limit
                self.line_limits[i] = max(self.line_limits[i], len(line))
                # pad the line with spaces if necessary
                lines[i] = f"{line: <{self.line_limits[i]}}"
        # If there are fewer lines than before, print additional blank lines
        if len(lines) < len(self.line_limits):
            for i in range(len(lines), len(self.line_limits)):
                lines.append(" " * self.line_limits[i])
        # Join the lines into a single string and print it
        print("\n".join(lines), flush=True)  # noqa: T201


class TimeDelta(timedelta):
    """
    Wrapper class for timedelta objects for better string formatting.
    """

    def __format__(self, format_spec: str) -> str:
        """Format a TimeDelta object according to the given format specifier. %m will be replaced by the elapsed
        minutes, %s will be replaced by the elapsed seconds as a floating point number with 1 decimal place and %S
        will be replaced by the seconds rounded to the nearest integer. If minutes are displayed, seconds will be shown
        with leading zeros and only up to 60s.

        :param format_spec: format specifier (e.g. "%mm%Ss")
        :return: formatted string
        """
        if format_spec == "":
            format_spec = "%mm%Ss"

        # format patterns for minutes, seconds (float) and rounded seconds respectively
        # a preceding '%' disables the substitution
        minutes_pattern = re.compile(r"(?<!%)%m")
        seconds_pattern = re.compile(r"(?<!%)%s")
        rounded_seconds_pattern = re.compile(r"(?<!%)%S")

        minutes_format = "{minutes}"

        seconds = self.total_seconds()
        if minutes_pattern.search(format_spec) is not None:
            # if minutes are present the seconds can only go up to 60, therefore we display a leading zero if necessary
            rounded_seconds_format = "{rounded_seconds:02}"
            seconds_format = "{seconds:04.1f}"
            minutes, seconds = divmod(seconds, 60)
            minutes = round(minutes)
        else:
            rounded_seconds_format = "{rounded_seconds}"
            seconds_format = "{seconds:.1f}"
            minutes = None

        rounded_seconds = round(seconds)

        format_spec = re.sub(minutes_pattern, minutes_format, format_spec)
        format_spec = re.sub(seconds_pattern, seconds_format, format_spec)
        format_spec = re.sub(rounded_seconds_pattern, rounded_seconds_format, format_spec)
        return format_spec.format(minutes=minutes, seconds=seconds, rounded_seconds=rounded_seconds)


def process_iteration(i: int, convergence: List[str], dt: TimeDelta) -> str:
    """Processes the available information of a finished scf iteration into a string.
    The first line contains the number of the finished iteration and the time of the iteration.
    Any subsequent lines contain information about the convergence process if present.

    :param i: number of the finished iteration
    :param convergence: information about the convergence status
    :param dt: elapsed time encode in a TimeDelta object
    :return: formatted string containing all the available information
    """
    if i > 0:
        return (f"Last finished iteration: {i: 3} ({dt:%ss})\n" + "\n".join(c.strip() for c in convergence)).strip()
    return ""


class TimeOut:
    """Class for tracking the timeout of a process."""

    def __init__(self, cmd: str, start_time: float, max_running_time: int):
        """Initialize a class

        :param cmd: The command being executed.
        :param start_time: The start time (in seconds) of the execution.
        :param max_running_time: The maximum allowed run time of the command.
        """
        self.cmd = cmd
        self.start_time = start_time
        self.max_running_time = max_running_time

    def is_expired(self) -> bool:
        """Checks if the process has exceeded its maximum run time.

        :return: True if the process has exceeded its maximum run time otherwise returns False
        """
        return perf_counter() - self.start_time > self.max_running_time

    def check(self):
        """Raises a TimeoutExpired exception if the command has exceeded its maximum run time."""
        if self.is_expired():
            raise subprocess.TimeoutExpired(self.cmd, self.max_running_time)


async def info_stream(
    file: Path, process: subprocess.Popen, queue: asyncio.Queue, sentinel: Sentinel, timeout: TimeOut
):
    """
    This function reads the INFO.OUT line by line and searches for the start of a new SCF cycle and convergence
    progress. When a new cycle starts it forwards this information to the queue.
    Raises an ExcitingException if the process ends before the file reading is complete and a subprocess.TimeoutExpired
    exception if the process is running longer than the given timeout.

    :param file: The file to be read.
    :param process: The process that is being monitored.
    :param queue: The queue in which the iteration number and convergence data are put.
    :param sentinel: A marker indicating the end of the data.
    :param timeout: TimeOut object to check if process should be killed early.
    """

    # Regular expression to match "scf iteration number :<number>"
    scf_iteration = re.compile(r"iteration number\s*:\s*(\d+)", re.IGNORECASE)

    # Regular expression to match lines containing convergence target
    convergence_target = re.compile(r"\s*\w\s+\(target\)")

    # Wait until the file exists
    while not file.exists():
        await asyncio.sleep(0.1)

    iteration = -1  # Initialize iteration number
    convergence = []  # Initialize convergence data

    # Open the file for reading
    with open(file) as f:
        while True:
            timeout.check()

            line = f.readline()  # Read a line from the file
            # Check if a new line was written
            if line:
                # Check if the line matches the scf_iteration pattern
                iteration_match = scf_iteration.search(line)
                if iteration_match is not None:
                    # If a previous iteration was found, put the iteration number and convergence data into the queue
                    if iteration > 0:
                        await queue.put((iteration, convergence))
                    # Update the iteration number and reset the convergence data
                    iteration = int(iteration_match.group(1))
                    convergence = []
                    continue

                # Check if the line matches the convergence_target pattern
                convergence_match = convergence_target.search(line)
                if convergence_match is not None:
                    # If a match is found, add the line to the convergence data
                    convergence.append(line)
                    continue

                # If the line contains "stopped", put the final iteration number and convergence data into the queue,
                # followed by the sentinel
                if "stopped" in line:
                    await queue.put((iteration, convergence))
                    await queue.put((iteration, sentinel))
                    return
            else:
                # If the process has ended unexpectedly, raise an exception
                if process.poll() is not None:
                    raise ExcitingException
                # If the process is still running, wait for a short time before trying to read the next line
                await asyncio.sleep(0.1)


async def stream_consumer(queue: asyncio.Queue, sentinel: Sentinel):
    """
    This functions reads and processes the information extracted from the INFO.OUT, which are put in the queue.
    It prints the currents process to the terminal.

    :param queue: The queue from which the iteration number and convergence data are taken.
    :param sentinel: A marker indicating the end of the process.
    """
    # Initialize a MultiLinePrinter instance for printing iteration summaries
    printer = MultiLinePrinter()

    # Initialize iteration variables
    i = 0
    convergence = []
    iteration_start = perf_counter()
    finished_iteration_summary = ""

    # Continue consuming items until the sentinel value is encountered
    while convergence is not sentinel:
        # Wait for a short duration to avoid busy waiting
        await asyncio.sleep(0.1)

        # Check if the queue is not empty
        if not queue.empty():
            # Get the finished iteration number and the convergence info from the queue
            i, convergence = await queue.get()

            # If the sentinel value is encountered, return from the function
            if convergence is sentinel:
                return

            # Measure the elapsed time of this iteration and start a new measurement
            iteration_end = perf_counter()
            dt = TimeDelta(seconds=iteration_end - iteration_start)
            finished_iteration_summary = process_iteration(i, convergence, dt)
            iteration_start = iteration_end

        # Get the summary of the current iteration and print it
        # Add iteration number and elapsed time if convergence is not the sentinel value
        if convergence is sentinel:
            running_iteration_summary = ""
        else:
            running_iteration_summary = f"Iteration {i + 1} ({TimeDelta(seconds=perf_counter() - iteration_start):%ss})"

        iteration_summary = "\n".join((finished_iteration_summary, running_iteration_summary)).strip()
        printer.print(iteration_summary)


class MonitoredGroundStateRunner(BinaryRunner):
    """Class to execute an exciting groundstate calculation and print the process of the calculation to stdout."""

    async def async_run(self) -> SubprocessRunResults:
        """Run exciting.

        Doing all the checks as the parent class BinaryRunner. Then executes exciting and monitors asynchronously the
        INFO.OUT file to get the current status of the groundstate calculation.

        Special handling is performed if exciting stops unexpectedly or the time limit is reached.

        :return: the run results with output and error message, runner code and run time
        """
        execution_list, my_env = self.get_execution_list_and_env()

        # Before we can start to track the changes to the info file, we first need to delete it,
        # otherwise the info reader might read old lines and/or be stuck at a position that doesn't exist
        # because exciting itself resets the info file.
        directory = Path(self.directory)
        info_file = directory / "INFO.OUT"
        if info_file.exists():
            info_file.unlink()

        time_start = perf_counter()
        process = subprocess.Popen(
            execution_list, env=my_env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=directory
        )

        print("Exciting started!")  # noqa: T201
        # init the async queue for progress information in the info file
        queue = asyncio.Queue()
        # init the sentinel value, which indicates the end of the process
        sentinel = Sentinel()

        crashed = False
        timeout = TimeOut(" ".join(execution_list), time_start, self.time_out)
        try:
            await asyncio.gather(
                info_stream(info_file, process, queue, sentinel, timeout), stream_consumer(queue, sentinel)
            )
        except ExcitingException:
            crashed = True
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = (output.decode("utf-8") for output in process.communicate())
            stderr = "Binary Runner: Job timed out. \n\n" + stderr
            # alias for the newline character, because backslashes are not allowed inside expression
            # parts of f-strings
            newline = "\n"
            # print stdout and stderr as bold and red text
            print(  # noqa: T201
                f"Exciting timed out!\n\033[31;1m{newline.join((stdout.strip(), stderr.strip())).strip()}\033[0m"
            )
            return SubprocessRunResults(stdout, stderr, RunnerCode.time_out, self.time_out)

        stdout, stderr = (output.decode("utf-8") for output in process.communicate())
        total_time = perf_counter() - time_start
        if not crashed:
            print(f"Exciting Finished! ({TimeDelta(seconds=total_time)})")  # noqa: T201
            warning_file = directory / "WARNINGS.OUT"
            if warning_file.exists():
                print("Exciting produced the following warnings:")  # noqa: T201
                print(f"\033[33;1m{warning_file.read_text('utf-8').strip()}\033[0m")  # noqa: T201
        else:
            newline = "\n"
            # print stdout and stderr as bold and red text
            print(  # noqa: T201
                f"Exciting crashed!\n\033[31;1m{newline.join((stdout.strip(), stderr.strip())).strip()}\033[0m"
            )
        return SubprocessRunResults(stdout, stderr, process.returncode, total_time)

    def run(self) -> SubprocessRunResults:
        """Runs the async run function.

        :return: the run results with output and error message, runner code and run time
        """
        return asyncio.run(self.async_run())
