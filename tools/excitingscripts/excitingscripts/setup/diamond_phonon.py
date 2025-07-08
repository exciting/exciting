"""Diamond phonon setup script."""

import json
import shutil
import sys
from argparse import ArgumentParser
from pathlib import Path
from typing import Callable, Union

import numpy as np
from excitingtools import ExcitingInputXML
from numpy.testing import assert_allclose
from numpy.typing import NDArray
if sys.version_info >= (3, 8):
    from typing import Protocol
else:
    from typing_extensions import Protocol


class SetupFunc(Protocol):
    """Protocol for diamond phonon setup function, providing a flexible way of performing type checking."""

    def __call__(
        self,
        maximum_displacement: float,
        displacement_points: int,
        work_directory: str,
        phonon_mode: Union[str, None] = None,
    ) -> None: ...


def point_specific_setup(
    get_new_positions: Callable[[float, NDArray[np.float64], Union[str, None]], NDArray[np.float64]],
    set_supercell: Union[Callable[[ExcitingInputXML], None], None] = None,
) -> SetupFunc:
    """Setup function for diamond phonon calculations.

    :param get_new_positions: function to get the new positions
    :param set_supercell: how to set a supercell, only for X phonons
    :return: setup function
    """

    def setup_diamond_phonon(
        maximum_displacement: float,
        displacement_points: int,
        work_directory: Union[str, Path],
        phonon_mode: Union[str, None] = None,
    ) -> None:
        """Set up a phonon calculation for diamond at the gamma or x point.

        Assumes to have an input file in the current directory which contains a diamond structure.

        :param maximum_displacement: Maximum displacement of the atom(s) in the unit cell
        :param displacement_points: The number of displacement points, should be in the range of 5 to 99 and odd
        :param work_directory: Working directory
        :param phonon_mode: which phonon mode to compute, should be TA, LA, TO or LO
         Only for X phonon calculations. Should be 'None' for gamma point. Legend:
         TA: transverse acoustic mode
         LA: longitudinal acoustic mode
         TO: transverse optical mode
         LO: longitudinal optical mode
        """
        input_file = Path("input.xml")
        if not input_file.exists():
            raise FileNotFoundError("Input file input.xml not found")

        assert (
            0 < maximum_displacement < 1
        ), f"Maximum displacement must be in range of 0 to 1, but is {maximum_displacement}"
        assert (
            4 < displacement_points <= 99
        ), f"Number of displacements must be in range of 5 to 99, but is {displacement_points}"
        assert displacement_points % 2 == 1, f"Number of displacements must be odd."

        input_obj = ExcitingInputXML.from_xml(input_file)
        structure = input_obj.structure
        assert len(structure.unique_species) == 1, "only a single species in the unit cell is supported"
        assert structure.positions == [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]], "given input file is not diamond structure"
        assert_allclose(
            structure.lattice,
            np.array([[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]),
            err_msg="basevectors in the original input file should be fcc lattice",
        )

        work_directory = Path(work_directory)
        if work_directory.exists():
            shutil.rmtree(work_directory)
        work_directory.mkdir()

        scale = getattr(structure.crystal_properties, "scale", 1)
        assert_allclose(
            np.array(getattr(structure.crystal_properties, "stretch", [1, 1, 1])),
            np.array([1, 1, 1]),
            err_msg="the original input file shouldn't contain a stretch value",
        )

        info = {
            "Equilibrium lattice parameter (alat) in a.u.": scale,
            "Maximum displacement (u,u,u); u in alat": maximum_displacement,
            "Number of displacements": displacement_points,
            "Volume of equilibrium unit cell in (a.u)^3": abs(np.linalg.det(structure.get_lattice())),
        }

        minimum_displacement = 1e-6
        displacements = [
            round(x, 15) if abs(x) > minimum_displacement else minimum_displacement
            for x in np.linspace(-maximum_displacement, maximum_displacement, displacement_points)
        ]

        if phonon_mode:
            phonon_mode = phonon_mode.upper()
            assert phonon_mode in ["TA", "LA", "TO", "LO"], "phonon mode not allowed"
            info["X-phonon-calculation mode"] = phonon_mode
            assert set_supercell is not None, "Need supercell method for X calculations"
            set_supercell(input_obj)
            displacements = displacements[int((displacement_points - 1) / 2):]

        equilibrium_positions = np.array(structure.positions)

        with open(work_directory / "INFO-diamond-phonon.json", "w") as f:
            json.dump(info, f, indent=4)

        for displacement in displacements:
            displacement_directory = work_directory / f"displ_{displacement}"
            displacement_directory.mkdir(exist_ok=True)

            structure.positions = list(get_new_positions(displacement, equilibrium_positions, phonon_mode))
            input_obj.write(displacement_directory / "input.xml")

    return setup_diamond_phonon


def arg_parser(point: str) -> ArgumentParser:
    """Get the arg parser for phonon setup scripts.

    :param point: the point, gamma or x
    :return: the argparser
    """
    parser = ArgumentParser(
        description=f"Script to set up a series of calculations with displacement for diamond "
        f"to compute phonons at the {point} point."
    )

    parser.add_argument(
        "maximum_displacement", type=float, help="Maximum displacement of the atom(s) in the unit cell."
    )

    parser.add_argument(
        "displacement_points", type=int, help="The number of displacement points in range of 5 to 99"
    )

    parser.add_argument(
        "--work-directory",
        "-w",
        type=str,
        default="workdir",
        dest="work_directory",
        help="Working directory",
    )

    return parser


if __name__ == "__main__":
    print("This is no script to call directly. Please call a point-specific version!")
