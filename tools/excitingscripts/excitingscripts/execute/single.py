"""Run a single **exciting** calculation.

Located at `excitingscripts/execute/single.py`.

Call as:

```bash
python3 -m excitingscripts.execute.single -r rundir
```
Where <code>rundir</code> is an optional parameter which specifies the running directory. If <code>rundir</code> is not specified, the calculation will run in the directory where the script is called.
"""

import os
import pathlib
from argparse import ArgumentParser

from excitingtools.runner.runner import BinaryRunner, RunnerCode


def run_exciting(root_directory: str=os.getcwd(), 
                 excitingroot: str=os.getenv("EXCITINGROOT"), 
                 filename: str="input.xml", 
                 timeout: int=3000) -> None:
    """Execute an exciting calculation in a given running directory.

    :param root_directory: Root directory.
    :param excitingroot: Environment variable string.
    :param filename: Name of the exciting input file
    :param timeout: Maximum runtime in seconds
    """
    if not excitingroot:
        raise ValueError(
            "EXCITINGROOT is not defined as an environment variable in the shell.\n"
            "If using bash please type: `export EXCITINGROOT=<path-to-exciting_smp>`")

    binary = pathlib.Path(excitingroot) / "install/bin/exciting_smp"
    # downward compatibility with old installations of exciting
    if not os.path.exists(binary):
        binary = pathlib.Path(excitingroot) / "bin/exciting_smp"
    n_threads = os.cpu_count()
    n_threads = 4 if n_threads is None else n_threads
    runner = BinaryRunner(binary, omp_num_threads=n_threads, time_out=timeout, directory=root_directory, args=[filename])
    result = runner.run()

    if result.return_code == RunnerCode.time_out:
        raise TimeoutError("exciting runtime exceeded.")

    if not result.success:
        print("Standard out:", result.stdout)
        print("Standard error:", result.stderr)
        raise RuntimeError("Running exciting failed")


def main() -> None:
    parser = ArgumentParser(description="""Execute a single exciting calculation in a given running directory.""")

    parser.add_argument("--root-directory", "-r",
                        default=[os.getcwd()],
                        nargs=1,
                        dest="root_directory",
                        help="root path for files that are created by this script")

    parser.add_argument("--input-file", "-f",
                        default=["input.xml"],
                        nargs=1,
                        dest="input_file",
                        help="name of the exciting input file")
        
    parser.add_argument("--timeout", "-t",
                        default=[1200],
                        nargs=1,
                        dest="timeout",
                        help="maximum runtime of calculation in seconds")

    args = parser.parse_args()

    run_exciting(args.root_directory[0], filename=args.input_file[0], timeout=args.timeout[0])


if __name__ == "__main__":
    main()
