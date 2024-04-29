import os
from argparse import ArgumentParser

from excitingtools.runner.monitored_runner import MonitoredGroundStateRunner


def main():

    parser = ArgumentParser()
    parser.add_argument("arg", default=None, nargs="?")
    parser.add_argument("--omp-threads", type=int, default=os.cpu_count(), help="Number of OMP threads.")
    parser.add_argument(
        "--timeout",
        type=int,
        default=31536000,  # default is equivalent to one year
        help="Number of seconds before a job is defined to have timed out.",
    )
    parser.add_argument("--mpi-procs", "-np", type=int, default=1, help="Number of MPI processes.")

    executable = parser.add_mutually_exclusive_group()
    binaries = {"smp": "smp", "serial": "serial", "mpi": "mpi", "mpismp": "mpismp", "debug": "debug_serial"}
    for binary in binaries:
        executable.add_argument(f"--{binary}", action="store_true")

    args = parser.parse_args()
    exciting_binary = "exciting"
    for binary in binaries:
        if getattr(args, binary) is True:
            exciting_binary = f"{exciting_binary}_{binaries[binary]}"
            break

    run_cmd = ""
    if "mpi" in exciting_binary:
        run_cmd = f"mpirun -np {args.mpi_procs}"

    runner = MonitoredGroundStateRunner(
        exciting_binary, run_cmd, args.omp_threads, args.timeout, args=None if args.arg is None else [args.arg]
    )
    runner.run()


if __name__ == "__main__":
    main()
