import os
from argparse import ArgumentParser
from pathlib import Path
from typing import Union

from excitingtools import ExcitingInputXML, parse
from scipy.constants import physical_constants


def convert_q_points_from_atomic_units_to_inverse_cm(run_dir: Union[str, Path]):
    """Display frequencies of the phonon modes from the exciting output file PHONON.OUT in inverse centimeters instead
    of atomic units.

    :param run_dir: Root directory.
    """
    hartree_2_inverse_cm = physical_constants["hartree-inverse meter relationship"][0] * 10**(-2)

    input_file = f"{run_dir}/input.xml"
    parsed_input = ExcitingInputXML.from_xml(input_file)

    cartesian = getattr(parsed_input.structure, "cartesian", False)

    if cartesian:
        print("Atomic positions are in cartesian coordinates")

    num_atoms = len(parsed_input.structure.positions)

    print("Total number of atoms    :", ("%5i"%num_atoms))

    if not Path(f"{run_dir}/PHONON.OUT").exists():
        raise FileNotFoundError("PHONON.OUT file not found")

    phonon_data = parse(f"{run_dir}/PHONON.OUT")

    num_q_points = len(phonon_data)
    print("Total number of q-points    :", ("%5i" % num_q_points))

    for q_point_index, q_point_data in phonon_data.items():
        print('\nq-point', ('%3i' % int(q_point_index)), ' is ', q_point_data['q_vector'])

        for mode in q_point_data['modes']:
            frequency_cm1 = mode['frequency'] * hartree_2_inverse_cm
            print('   mode', ('%3i' % int(mode['mode_index'])), '     frequency:', '%10.4f cm-1' % frequency_cm1)


def main() -> None:
    parser = ArgumentParser(description="""Display frequencies of the phonon modes from the exciting output file
                                        PHONON.OUT in inverse centimeters instead of atomic units.""")

    parser.add_argument("--root-directory", "-r",
                        default=[os.getcwd()],
                        nargs=1,
                        dest="root_directory",
                        help="root path containing files needed to run this script")

    args = parser.parse_args()

    convert_q_points_from_atomic_units_to_inverse_cm(args.root_directory[0])


if __name__ == "__main__":
    main()
