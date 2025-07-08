"""Python script for generating strained structures."""

import os
import pathlib
import numpy as np
from argparse import ArgumentParser
from os.path import join
from typing import Union
from excitingtools import ExcitingInputXML
import shutil
import json


def setup_displaced_structures(input_file: Union[str, pathlib.Path], eta_strain: float, umax: float,
                               number_of_displ: int, workdir=os.getcwd()) -> None:
    """
    :param input_file: Input file.
    :param eta_strain: Lagrangian strain.
    :param umax: Maximum displacement.
    :param number_of_displ: Number of displacements.
    :param workdir: Working directory.
    """
    parsed_input = ExcitingInputXML.from_xml(input_file)

    # Extract crystal properties
    if hasattr(parsed_input.structure.crystal_properties, "scale"):
        ref_scale = parsed_input.structure.crystal_properties.scale
    else:
        raise ValueError("\nERROR: There is NO scale attribute in input.xml!")

    if hasattr(parsed_input.structure.crystal_properties, "stretch"):
        stretch = parsed_input.structure.crystal_properties.stretch
    else:
        stretch = [1.0, 1.0, 1.0]

    base_vectors = np.array(parsed_input.structure.lattice)

    # Create work directory
    os.makedirs(workdir, exist_ok=True)

    if umax == 0:
        parsed_input.write(join(workdir, "input.xml"))

        # Write strain value to file
        with open(join(workdir, "displ"), "w") as output_str:
            output_str.write("0.00")

        print("Single un-displaced calculation\n")
        return

    # Calculate initial volume
    volume = abs(np.linalg.det(base_vectors) * ref_scale ** 3 * np.prod(stretch))

    eps_strain = -1.0 + np.sqrt(1.0 + 2.0 * eta_strain)
    new_scale = (1.0 + eps_strain) * ref_scale

    with open(join(workdir, 'INFO-graphene-along-c'), "w") as output_info:
        output_info.write(
            f"\nLattice parameter (alat)              = {new_scale:11.8f} [a.u.]\n"
            f"Maximum displacement (0,0,umax); umax =  {umax} [c/a]\n"
            f"Number of displacements               =  {number_of_displ}\n"
            f"Lagrangian strain                     =  {eta_strain:11.8f}\n"
        )

    info = {
        "Equilibrium lattice parameter (alat) in a.u.": ref_scale,
        "Maximum displacement (u,u,u); u in alat": umax,
        "Number of displacements": number_of_displ,
        "Volume of equilibrium unit cell in (a.u)^3": volume,
    }

    phonon_mode = None
    info["X-phonon-calculation mode"] = phonon_mode

    # Parse the atomic positions
    atomic_positions = np.array(parsed_input.structure.positions)

    # Create a json file
    info_file_json = join(workdir, "INFO-diamond-phonon.json")
    with open(info_file_json, 'w') as f:
        json.dump(info, f, indent=4)

    # Creating source.xml file
    shutil.copy(input_file, join(workdir, 'source.xml'))

    # Info file
    with open(join(workdir, 'INFO-graphene-along-c'), "w") as output_info:
        output_info.write(
            f"\nLattice parameter (alat)              = {new_scale:11.8f} [a.u.]\n"
            f"Maximum displacement (0,0,umax); umax =  {umax} [c/a]\n"
            f"Number of displacements               =  {number_of_displ}\n"
            f"Lagrangian strain                     =  {eta_strain:11.8f}\n"
        )

    # Calculate strain steps
    delta = number_of_displ - 1

    if number_of_displ <= 1:
        delta = 1

    displ_step = umax / delta

    for i in range(number_of_displ):
        # Update displ value
        displ = i * displ_step

        # Tolerance on displ values
        if abs(displ) < 0.000001:
            displ = 0.000001
        if abs(displ - 0.50) < 0.000001:
            displ = 0.499999

        # Displacement matrix
        displ_matrix = np.zeros(shape=np.shape(atomic_positions))
        displ_matrix[-1][-1] = displ_matrix[-1][-1] + displ

        # Update the base vectors
        new_atomic_positions = atomic_positions.copy() + displ_matrix

        # Update base vectors in the XML
        parsed_input.structure.positions = new_atomic_positions.tolist()

        # Create rundir
        rundir = join(workdir, f"displ-{displ}")
        os.makedirs(rundir, exist_ok=True)

        # Write deformed structure to file
        parsed_input.write(join(rundir, "input.xml"))



def main() -> None:
    parser = ArgumentParser(description="Python script for generating strained structures.")

    parser.add_argument("--input-file", "-i",
                        type=str,
                        default=["input.xml"],
                        nargs=1,
                        dest="infile",
                        help="name of the input file")

    parser.add_argument("--work-directory", "-w",
                        type=str,
                        default=[join(os.getcwd(), "workdir")],
                        nargs=1,
                        dest="workdir",
                        help="path for folders that are created by this script")

    parser.add_argument("--lstrain", '--ls',
                        type=float,
                        required=True,
                        nargs=1,
                        dest="eta_strain",
                        help="The langrangian strain in the range [-0.5,0.5]")

    parser.add_argument("--maximum-displacement", '--u-max',
                        type=float,
                        required=True,
                        nargs=1,
                        dest="umax",
                        help="Maximum displacement(umax) in the range [0-1]")

    parser.add_argument("--number-of-displacements", '-n',
                        type=int,
                        required=True,
                        nargs=1,
                        dest="number_of_displ",
                        help="The number of displacements in [0,umax] having the range [1-99]")

    args = parser.parse_args()

    infile = join(os.getcwd(), args.infile[0])
    workdir = args.workdir[0]
    eta_strain = args.eta_strain[0]
    umax = args.umax[0]
    number_of_displ = args.number_of_displ[0]

    # Check for input file existence
    if not os.path.exists(infile):
        raise FileNotFoundError(infile)

    # Check lagrangian strain is in correct range
    if not abs(eta_strain) <= 0.5:
        print("\nERROR: Lagrangian strain is out of range [-0.5,0.5]!")
        return

    # Check umax is in correct range
    if not (0 <= umax <= 1):
        print("\nERROR: Maximum displacement is out of range [0-1]!")
        return

    # Check number of displacements is in correct range
    if not (1 <= number_of_displ <= 99):
        print("\nERROR: Number of displacements is out of range [1-99]!")
        return

    # Set up the deformed structures
    setup_displaced_structures(infile, eta_strain, umax, number_of_displ, workdir)



if __name__ == "__main__":
    main()
