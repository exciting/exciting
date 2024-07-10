"""Python visualization tool for following relative atomic coordinates of atoms during the relaxation process."""

from os.path import join
import math
import numpy as np
from argparse import ArgumentParser
from typing import List, Dict
from excitingtools import parse
from excitingscripts.utils.utils import get_structure_optimizations_properties, get_num_atoms, is_coordinate_cartesian
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ptk
import warnings

if matplotlib.__version__.split(".")[0] == "2":
    matplotlib.style.use('classic')


def calculate_bond_lengths(position_data: List[Dict], atom1: int, atom2: int, lattice_matrix: np.ndarray,
                           isCartesian=False, threshold=0.67) -> \
        List[float]:
    """ Calculate bond lengths between specified atoms from position data.

    :param position_data: list of dictionaries containing atomic positions.
    :param atom1: first atom number.
    :param atom2: second atom number.
    :param lattice_matrix: unit cell containing the molecule
    :param isCartesian: coordinate type, either "lattice" - False or "cartesian" - True.
    :param threshold: threshold value for periodic boundary conditions

    :return: list of calculated bond lengths.
    """
    bond_lengths = []
    delta_x, delta_y, delta_z = [], [], []
    iter = 0

    # Calculation of cell dim (considering orthogonal)
    cell_dim = [1.0, 1.0, 1.0]

    # The primitive lattice vectors should be along the x,y,z axis
    if isCartesian:
        if lattice_matrix is None:
            warnings.warn("No Lattice matrix found despite cartesian system")
        else:
            diagonals = np.diagonal(lattice_matrix)
            if np.any(np.abs(lattice_matrix - np.diag(diagonals))) > 1e-8:
                raise ValueError("Lattice type non implemented!")
            cell_dim = diagonals.tolist()

    for positions in position_data:
        pos1 = positions[atom1]
        pos2 = positions[atom2]

        # If not the first iteration, check for periodic boundary conditions
        delta_r = [pos2[i] - pos1[i] for i in range(3)]
        if iter > 0:
            prev_delta_r = [delta_x[iter - 1], delta_y[iter - 1], delta_z[iter - 1]]
            for i in range(3):
                if (delta_r[i] - prev_delta_r[i]) > threshold * cell_dim[i]:
                    delta_r[i] = delta_r[i] - cell_dim[i]
                if (prev_delta_r[i] - delta_r[i]) > threshold * cell_dim[i]:
                    delta_r[i] = delta_r[i] + cell_dim[i]

        # Append the adjusted differences to the lists
        delta_x.append(delta_r[0])
        delta_y.append(delta_r[1])
        delta_z.append(delta_r[2])

        # Increment of the iteration counter
        iter += 1

        # Find the maximum deviation component
        max_abs_delta = abs(max(delta_r, key=abs))

        # Normalize data
        initial_delta = delta_x[0]
        if max_abs_delta < 1000:
            initial_delta = 0.0

        for i in range(len(delta_x)):
            delta_x[i] -= initial_delta
            delta_y[i] -= initial_delta
            delta_z[i] -= initial_delta

        # Adjust the values that exceed half the cell dimensions
        for i in range(len(delta_x)):
            if delta_x[i] >= cell_dim[0] / 2.:
                delta_x[i] = cell_dim[0] - delta_x[i]
            if delta_y[i] >= cell_dim[1] / 2.:
                delta_y[i] = cell_dim[1] - delta_y[i]
            if delta_z[i] >= cell_dim[2] / 2.:
                delta_z[i] = cell_dim[2] - delta_z[i]

    for i in range(len(delta_x)):
        bond_length = math.sqrt(abs(delta_x[i] ** 2 + delta_y[i] ** 2 + delta_z[i] ** 2))
        bond_lengths.append(bond_length)

    return bond_lengths


def plot_bond_lengths(bond_lengths: List[float], run_dir: str, isCartesian: bool, atom1: int, atom2: int, show: bool,
                      dpi: int):
    """ Plot the bond lengths and save the plot to a file.

    :param bond_lengths: list of bond lengths to plot.
    :param run_dir: directory where exciting runs.
    :param isCartesian: coordinate type, either "lattice" - False or "cartesian" - True.
    :param atom1:  first atom number
    :param atom2: second atom number
    :param show: whether to display the plot.
    :param dpi: resolution in DPI for the saved plot.
    """
    plt.rcParams.update({
        'figure.figsize': (10, 7.5),
        'ytick.minor.size': 6,
        'xtick.major.pad': 8,
        'ytick.major.pad': 4,
        'patch.linewidth': 2.0,
        'axes.linewidth': 4.0,
        'lines.markersize': 10,
        'lines.linewidth': 2.0,
        'axes.formatter.limits': (-4, 4),
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'axes.labelsize': 20,
        'legend.fontsize': 25,
        'mathtext.default': 'regular'})

    fig, ax = plt.subplots()

    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(2)

    acoord = "cartesian" if isCartesian else "lattice"

    ax.plot(bond_lengths, 'ro-', label=f'Bond {atom1}-{atom2}', )
    ax.set_xlabel('Optimization steps')
    ax.set_ylabel(f'Bond length ({acoord}) [bohr]')
    ax.grid(True, linestyle='--')

    ax.xaxis.set_major_locator(ptk.MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(ptk.ScalarFormatter(useOffset=True, useMathText=True))

    # Calculate x-axis limits
    iter_count = len(bond_lengths)
    xmin = 0 - (iter_count - 1) / 20.0
    xmax = (iter_count - 1) + (iter_count - 1) / 20.0

    if iter_count == 1:
        xmin = -1
        xmax = 1

    # Set x-axis limits
    ax.set_xlim(xmin, xmax)

    # Adjust legend
    plt.legend(borderaxespad=.8)

    plt.savefig(join(run_dir, "PLOT.png"), format='png', dpi=dpi)

    if show:
        plt.show()

    plt.close()


def main() -> None:
    parser = ArgumentParser(
        description="Python visualization tool for following relative atomic coordinates of atoms during the relaxation process.")

    parser.add_argument("directory",
                        type=str,
                        help="Directory where exciting runs")

    parser.add_argument("atom1",
                        type=int,
                        help="First atom number")

    parser.add_argument("atom2",
                        type=int,
                        help="Second atom number")

    parser.add_argument("-sh", "--show",
                        action="store_true",
                        help="Show plot in a window")

    parser.add_argument("--dpi",
                        type=int,
                        default=300,
                        help="Resolution in DPI for saved plot")

    args = parser.parse_args()
    run_dir = args.directory
    atom1 = args.atom1
    atom2 = args.atom2

    # Get the coordinate type
    isCartesian = is_coordinate_cartesian(run_dir)

    # Get the number of atoms
    num_atoms = get_num_atoms(run_dir)

    # Check if the specified atom numbers exceed the total number of atoms
    if atom1 > num_atoms or atom2 > num_atoms:
        raise ValueError(f"The specified atoms {atom1} or {atom2} exceed the total number of atoms {num_atoms}.")

    # Find the lattice matrix
    lattice_matrix = None
    try:
        input_file = join(run_dir, "input.xml")
        parsed_input = parse(input_file)
        lattice_matrix = np.array(parsed_input['structure']['lattice'])
    except KeyError:
        pass

    # Extract data from the file
    position_data: List[Dict] = get_structure_optimizations_properties(run_dir, "Atomic positions")

    # Calculate bond lengths
    bond_lengths = calculate_bond_lengths(position_data, atom1, atom2, lattice_matrix, isCartesian)

    # Plot bond lengths
    plot_bond_lengths(bond_lengths, run_dir, isCartesian, atom1, atom2, args.show, args.dpi)


if __name__ == "__main__":
    main()
