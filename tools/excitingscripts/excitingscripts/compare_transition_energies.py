import pathlib
from argparse import ArgumentParser
from typing import Union, Tuple

from excitingtools.exciting_obj_parsers.eigenvalue_parser import parse_eigenvalues
from scipy.constants import physical_constants


def determine_transition_energies(root_directory: Union[str, pathlib.Path]) -> Tuple[float, float]:
    """Determine the transition energies for the transitions Γ→Γ and Γ→X for a given directory in which an exciting
    calculation was performed.

    :param root_directory: Root directory.
    """

    ha_to_ev = physical_constants["hartree-electron volt relationship"][0]
    eigval_data= parse_eigenvalues(f"{root_directory}/eigval.xml")
    eigval_data.all_eigenvalues = eigval_data.all_eigenvalues * ha_to_ev

    gamma = [0, 0, 0]
    X = [0.5, 0.5, 0]

    gap_gamma_gamma = eigval_data.get_transition_energy(gamma, gamma)
    gap_gamma_X = eigval_data.get_transition_energy(gamma, X)

    return gap_gamma_gamma, gap_gamma_X

def main() -> None:
    parser = ArgumentParser(description="Determine the transition energies for the transitions Γ→Γ and Γ→X.")


    parser.add_argument("--root-directories", "-r",
                        nargs='+',
                        dest="root_directories",
                        help="names of root directories")

    args = parser.parse_args()

    transition_energies = [[""], ["Gamma -> Gamma:"],["Gamma -> X:"]]
    transition_energies[0].extend(args.root_directories)
    for i in range(len(args.root_directories)):
        transition_energies[1].append(str(round(determine_transition_energies(args.root_directories[i])[0], 3)))
        transition_energies[2].append(str(round(determine_transition_energies(args.root_directories[i])[1], 3)))

    print("\n------------------------------------------------\n")
    print(" Transition energies in eV:\n")
    n = max(len(x) for entry in transition_energies for x in entry)
    for entry in transition_energies:
        print(''.join(x.ljust(n + 2) for x in entry))
    print("\n------------------------------------------------\n")


if __name__ == "__main__":
    main()
