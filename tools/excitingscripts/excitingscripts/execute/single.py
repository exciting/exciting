"""Run a single **exciting** calculation.

Located at `excitingscripts/execute/single.py`.

Call as:

```bash
python3 -m excitingscripts.execute.single -r rundir
```
Where <code>rundir</code> is an optional parameter which specifies the running directory. If <code>rundir</code> is not specified, the calculation will run in the directory where the script is called.
"""

import os
import sys
import re
import pathlib
import psutil
from argparse import ArgumentParser

from excitingtools.runner.runner import BinaryRunner, RunnerCode


def run_exciting(root_directory: str=os.getcwd(), 
                 excitingroot: str=os.getenv("EXCITINGROOT"), 
                 filename: str="input.xml", 
                 timeout: int=3000,
                 verbose: bool = False, 
                 verify: bool = False) -> None:
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

    binary_path = None
    for binary_name in ['exciting_smp', 'exciting_mpismp']:
        for install_dir in ['install/bin', 'bin']:
            binary_path = pathlib.Path(excitingroot) / f"{install_dir}/{binary_name}"
            if not os.path.exists(binary_path):
                binary_path = None
            else:
                break
        if binary_path is not None:
            break
    n_threads = psutil.cpu_count(logical=False)
    n_threads = 4 if n_threads is None else n_threads
    n_threads = min([n_threads, 8])
    runner = BinaryRunner(binary_path, omp_num_threads=n_threads, time_out=timeout, directory=root_directory, args=[filename])
    result = runner.run()

    stopped_INFO_OUT = True

    if result.return_code == RunnerCode.time_out:
        raise TimeoutError("exciting runtime exceeded.")

    if verify:
        with open(os.path.join(root_directory, 'INFO.OUT'), 'r') as f_:
            info_out = f_.read()
        if re.search(r'EXCITING\s+[A-Z]+\sstopped', info_out) is None:
            stopped_INFO_OUT = False

    if verbose or not result.success or not stopped_INFO_OUT:
        if len(result.stdout) != 0:
            print("STDOUT:\n", result.stdout, file=sys.stdout)
        if len(result.stderr) != 0:
            print("STDERR:\n", result.stderr, file=sys.stderr)

    assert stopped_INFO_OUT, "INFO.OUT does not include exciting stop message. Likely an error has occured."

    if not result.success:
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
                        default=[3000],
                        nargs=1,
                        type=int,
                        dest="timeout",
                        help="maximum runtime of calculation in seconds")

    parser.add_argument("--verbose", "-v",
                        dest="verbose",
                        action="store_true",
                        help="print STDOUT and STDERR if present")

    parser.add_argument("--verify",
                        dest="verify",
                        action="store_true",
                        help="check for the message that exciting stopped in INFO.OUT")

    args = parser.parse_args()

    run_exciting(args.root_directory[0], 
                 filename=args.input_file[0], 
                 timeout=args.timeout[0], 
                 verbose=args.verbose, 
                 verify=args.verify)


if __name__ == "__main__":
    main()
