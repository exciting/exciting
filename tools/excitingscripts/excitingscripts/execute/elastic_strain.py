"""Run a series of **exciting** calculations with different strain values.

Located at `excitingscripts/execute/elastic_strain.py`.

Call as:

```bash
python3 -m excitingscripts.execute.elastic_strain 
```
"""

import os
import pathlib
from argparse import ArgumentParser
from typing import Union

import numpy as np
import numpy.typing as npt
from excitingscripts.execute.single import run_exciting
from excitingtools import parse
from scipy.constants import physical_constants


def execute_elastic_strain(
    root_directory: Union[str, pathlib.Path] = os.getcwd(),
    dft_half: bool = False,
    excitingroot=os.getenv("EXCITINGROOT"),
) -> npt.NDArray[np.float64]:
    """Execute a series of exciting calculations with different interlayer distances.

    :param root_directory: Root directory.
    :param dft_half: Boolean with "True" value for DFT-1/2 calculations.
    :param excitingroot: Environment variable string.
    :returns: Array with energy-strain data.
    """

    # Counting the number of rundir-xx
    listdir = os.listdir(root_directory)
    displ_points = len(listdir)

    for directory in listdir:
        if not directory.startswith("rundir"):
            displ_points -= 1

    rundir_infty = f"{root_directory}/rundir-oo"
    if os.path.exists(rundir_infty):
        displ_points -= 1

    ha_to_ev = physical_constants["hartree-electron volt relationship"][0]

    energy = []
    strain_values = []

    for i in range(displ_points):
        run_exciting(f"{root_directory}/rundir-{i + 1}", excitingroot)
        results = parse(f"{root_directory}/rundir-{i + 1}/INFO.OUT")
        max_scf = max([int(j) for j in results["scl"].keys()])
        converged_results = results["scl"][str(max_scf)]
        if dft_half:
            energy.append(converged_results["Estimated fundamental gap"] * ha_to_ev)
        else:
            energy.append(converged_results["Total energy"])

        with open(f"{root_directory}/rundir-{i + 1}/strain-{i + 1}") as f:
            strain_values.append(float(f.readline()))

    # Only valid for "normal" elastic-strain calculations, not DFT-1/2
    if os.path.exists(rundir_infty):
        run_exciting(rundir_infty, excitingroot, "input.xml", 1500)

        results = parse(rundir_infty + "/INFO.OUT")
        max_scf = max([int(i) for i in results["scl"].keys()])
        converged_results = results["scl"][str(max_scf)]
        energy.append(converged_results["Total energy"])

        with open(rundir_infty + "/strain-oo") as f:
            strain_values.append(float(f.readline()))

    energy_strain_data = np.vstack((strain_values, energy)).T

    return energy_strain_data


def main() -> None:
    parser = ArgumentParser(
        description="""Execute a series of exciting calculations."""
    )

    parser.add_argument(
        "--root-directory",
        "-r",
        default=[os.getcwd()],
        nargs=1,
        dest="root_directory",
        help="root path for files that are created by this script",
    )

    parser.add_argument(
        "--dft-half",
        dest="dft_half",
        action="store_true",
        help=""" If present, data will be saved to "bandgap-vs-rcut".
                             Otherwise, data will be saved to "energy-vs-strain".""",
    )

    parser.set_defaults(dft_half=False)

    args = parser.parse_args()

    energy_strain_data = execute_elastic_strain(args.root_directory[0], args.dft_half)

    output_file = "bandgap-vs-rcut" if args.dft_half else "energy-vs-strain"

    with open(f"{args.root_directory[0]}/{output_file}", "w") as f:
        np.savetxt(f, energy_strain_data, fmt="%15.8f %15.10f")


if __name__ == "__main__":
    main()
