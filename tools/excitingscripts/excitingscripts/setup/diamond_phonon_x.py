"""Script to set up a phonon calculation for diamond at the X point."""

import numpy as np
from excitingtools import ExcitingInputXML
from numpy.typing import NDArray

from excitingscripts.setup.diamond_phonon import arg_parser, point_specific_setup


def get_new_positions(d: float, equilibrium_positions: NDArray[np.float64], phonon_mode: str) -> NDArray[np.float64]:
    """Get the new displaced positions for a phonon with X character.

    :param d: displacement of the atoms
    :param equilibrium_positions: the non-displaced positions
    :param phonon_mode: Mode of the phonon
    :return: array with the new positions
    """
    z = np.sqrt(2)
    mode_displacement_map = {
        "LA": [[0.0, 0.0, d], [0.0, 0.0, d], [0.0, 0.0, -d], [0.0, 0.0, -d]],
        "LO": [[0.0, 0.0, d], [0.0, 0.0, -d], [0.0, 0.0, -d], [0.0, 0.0, d]],
        "TA": [[0.0, z * d, 0.0], [0.0, z * d, 0.0], [0.0, -z * d, 0.0], [0.0, -z * d, 0.0]],
        "TO": [[0.0, z * d, 0.0], [0.0, -z * d, 0.0], [0.0, -z * d, 0.0], [0.0, z * d, 0.0]],
    }
    return equilibrium_positions + np.array(mode_displacement_map[phonon_mode])


def set_supercell(input_obj: ExcitingInputXML) -> None:
    """Set a supercell according to an x phonon.

    :param input_obj: the input xml object
    """
    structure = input_obj.structure

    # add two more atoms and create a supercell
    species = structure.unique_species[0]
    for _ in range(2):
        structure.add_atom(species, [0, 0, 0])
    structure.positions = [[0.0, 0.0, 0.0], [0.0, 0.5, 0.25], [0.5, 0.5, 0.5], [0.5, 1.0, 0.75]]

    scale = getattr(structure.crystal_properties, "scale", 1)
    structure.crystal_properties.scale = scale / np.sqrt(2)
    structure.lattice = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, np.sqrt(2)]])

    factors = [1, 1, 1 / np.sqrt(2)]
    new_ngridk = [int(x * factors[i] / 2 ** (1 / 6) + 0.5) for i, x in enumerate(input_obj.groundstate.ngridk)]
    input_obj.groundstate.ngridk = new_ngridk


def main() -> None:
    parser = arg_parser("X")

    parser.add_argument(
        "--phonon-mode",
        "-p",
        type=str,
        dest="phonon_mode",
        help="Phonon mode to compute, either TA, LA, TO or LO",
        required=True,
    )

    args = parser.parse_args()

    setup_func = point_specific_setup(get_new_positions, set_supercell)
    setup_func(args.maximum_displacement, args.displacement_points, args.work_directory, args.phonon_mode)


if __name__ == "__main__":
    main()
