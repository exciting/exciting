import os
import pathlib
from argparse import ArgumentParser
from os.path import join
from typing import Union

import numpy as np
from excitingscripts.execute.single import run_exciting

from excitingtools import parse


def execute_convergence_test(k_initial: int, k_final: int, rgkmax_initial: int,
                             rgkmax_final: int, root_directory=os.getcwd(),
                             excitingroot=os.getenv("EXCITINGROOT")) -> list:
    """Execute a series of exciting calculations with varying values for the groundstate attributes ngridk and rgkmax
     and return a list containing the total energy value for each set of parameters.

    :param k_initial: Initial k-value for defining the groundstate attribute ngridk.
    :param k_final: Final k-value for defining the groundstate attribute ngridk.
    :param rgkmax_initial: Initial value for the groundstate attribute rgkmax.
    :param rgkmax_final: Final value for the groundstate attribute rgkmax.
    :param root_directory: Root directory.
    :param excitingroot: Environment variable string.
    :returns: List containing total energy values for each set of parameters.
    """

    dk = 2
    k_final = k_final + 1
    rgkmax_final = rgkmax_final + 1

    total_energy_ngridk_rgkmax = np.empty([int((k_final - k_initial + 1) / dk), int(rgkmax_final - rgkmax_initial)])
    convergence_test = []

    for i_k, k in enumerate(range(k_initial, k_final, dk)):
        for i_rgkmax, rgkmax in enumerate(range(rgkmax_initial, rgkmax_final)):

            run_exciting(f"{root_directory}/{k}_{rgkmax}", excitingroot)

            results = parse(join(os.getcwd(), f"{root_directory}/{k}_{rgkmax}/INFO.OUT"))
            max_scf = max([int(i) for i in results["scl"].keys()])
            converged_results = results["scl"][str(max_scf)]
            total_energy_ngridk_rgkmax[i_k, i_rgkmax] = converged_results["Total energy"]
            convergence_test.append([k, rgkmax, total_energy_ngridk_rgkmax[i_k, i_rgkmax]])

    return convergence_test


def main() -> None:
    parser = ArgumentParser(description="""Execute a series of exciting calculations with varying values for the
                                        groundstate attributes ngridk and rgkmax.""")

    parser.add_argument("k_initial",
                        type=int,
                        nargs=1,
                        help="initial k-value for defining ngridk")

    parser.add_argument("k_final",
                        type=int,
                        nargs=1,
                        help="final k-value for defining ngridk")

    parser.add_argument("rgkmax_initial",
                        type=int,
                        nargs=1,
                        help="initial k-value for rgkmax")

    parser.add_argument("rgkmax_final",
                        type=int,
                        nargs=1,
                        help="final k-value for rgkmax")

    parser.add_argument("--root-directory", "-r",
                        type=Union[str, pathlib.Path],
                        default=[os.getcwd()],
                        nargs=1,
                        dest="root_directory",
                        help="root path for files that are created by this script")

    args = parser.parse_args()

    convergence_test = execute_convergence_test(args.k_initial[0], args.k_final[0], args.rgkmax_initial[0],
                                                args.rgkmax_final[0], args.root_directory[0])

    with open(f"{args.root_directory[0]}/convergence-test", "w") as f:
        np.savetxt(f, convergence_test, fmt='%4i  %6.2f  %18.8f')


if __name__ == "__main__":
    main()
