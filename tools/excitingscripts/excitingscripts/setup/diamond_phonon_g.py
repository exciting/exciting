"""Script to set up a phonon calculation for diamond for a phonon at the gamma point."""

from typing import Union

import numpy as np
from numpy.typing import NDArray

from excitingscripts.setup.diamond_phonon import arg_parser, point_specific_setup


def get_new_positions(
        displacement: float, equilibrium_positions: NDArray[np.float64], phonon_mode: Union[str, None]
) -> NDArray[np.float64]:
    """Get the new displaced positions for a phonon with Gamma character.

    :param displacement: displacement of the atom(s)
    :param equilibrium_positions: the non-displaced positions
    :param phonon_mode: not used here
    :return: array with the new positions
    """
    assert phonon_mode is None, "phonon mode not used for the gamma point."
    displacement_array = np.array([[0.0] * 3, [displacement] * 3])
    return equilibrium_positions + displacement_array


def main() -> None:
    args = arg_parser("Gamma").parse_args()

    setup_func = point_specific_setup(get_new_positions)
    setup_func(args.maximum_displacement, args.displacement_points, args.work_directory)


if __name__ == "__main__":
    main()
