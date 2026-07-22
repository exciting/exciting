"""Run a series of **exciting** calculations with different values of the main computational parameters.

Located at `excitingscripts/execute/convergence_test.py`.

Call as:

```bash
python3 -m excitingscripts.execute.convergence_test k_i k_f rgkmax_i rgkmax_f
```
Where <code>k_i</code> and <code>k_f</code> are the initial and final k-values for defining the <code><span style="color:green">groundstate</span></code> attribute <code><span style="color:MediumBlue">ngridk</span></code>, and <code>rgkmax_i</code> and <code>rgkmax_f</code> the initial and final values for the <code><span style="color:green">groundstate</span></code> attribute <code><span style="color:MediumBlue">rgkmax</span></code>.
"""

import os
from argparse import ArgumentParser
from os.path import join

import numpy as np
from excitingscripts.execute.single import run_exciting

from excitingtools import parse


def execute_convergence_test(k_initial: int, 
                             k_final: int,
                             delta_k: int, 
                             rgkmax_initial: int,
                             rgkmax_final: int,
                             delta_rgkmax: int, 
                             root_directory=os.getcwd(),
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

    k_final = k_final + 1
    rgkmax_final = rgkmax_final + (1 * delta_rgkmax)

    #total_energy_ngridk_rgkmax = np.empty([int((k_final - k_initial + 1) / delta_k), int(rgkmax_final - rgkmax_initial)])
    convergence_test = []

    for k in range(k_initial, k_final, delta_k):
        for rgkmax in np.arange(rgkmax_initial, rgkmax_final, delta_rgkmax):

            #print(f"Running for rgkmax {rgkmax} and ngridk {k}")
            run_exciting(f"{root_directory}/{k}_{rgkmax}", excitingroot)

            results = parse(join(os.getcwd(), f"{root_directory}/{k}_{rgkmax}/INFO.OUT"))
            max_scf = max([int(i) for i in results["scl"].keys()])
            converged_results = results["scl"][str(max_scf)]
            convergence_test.append([k, rgkmax, converged_results["Total energy"]])

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
    
    parser.add_argument('--delta-k', "-dk",
                        type=int,
                        default=2,
                        
                        help="change of k-grid value")

    parser.add_argument("rgkmax_initial",
                        type=int,
                        nargs=1,
                        help="initial k-value for rgkmax")

    parser.add_argument("rgkmax_final",
                        type=int,
                        nargs=1,
                        help="final k-value for rgkmax")
    
    parser.add_argument('--delta-rgkmax', "-dr",
                        type=float,
                        default=1,
                        help="change of k-grid value")

    parser.add_argument("--root-directory", "-r",
                        type=str,
                        default=[os.getcwd()],
                        nargs=1,
                        dest="root_directory",
                        help="root path for files that are created by this script")

    args = parser.parse_args()

    convergence_test = execute_convergence_test(args.k_initial[0], 
                                                args.k_final[0], 
                                                args.delta_k, 
                                                args.rgkmax_initial[0],
                                                args.rgkmax_final[0], 
                                                args.delta_rgkmax, 
                                                args.root_directory[0])

    with open(f"{args.root_directory[0]}/convergence-test", "w") as f:
        np.savetxt(f, convergence_test, fmt='%4i  %6.2f  %18.8f')


if __name__ == "__main__":
    main()
