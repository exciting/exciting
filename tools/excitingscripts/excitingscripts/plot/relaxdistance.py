"""Python visualization tool for following relative atomic coordinates of atoms during the relaxation process."""
import os
import warnings
from argparse import ArgumentParser
from os.path import join
from typing import List, Dict

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ptk
import numpy as np
from excitingscripts.utils.utils import get_structure_optimizations_properties, get_num_atoms, is_coordinate_cartesian
from excitingtools import parse

if matplotlib.__version__.split(".")[0] == "2":
    matplotlib.style.use('classic')


def calculate_relative_coordinates(position_data: List[Dict], atom1: int, atom2: int, lattice_matrix: np.ndarray,
                                   isCartesian=False, threshold=0.67) -> List[float]:
    """ Calculate bond lengths between specified atoms from position data.

    :param position_data: list of dictionaries containing atomic positions.
    :param atom1: first atom number.
    :param atom2: second atom number.
    :param lattice_matrix: unit cell containing the molecule
    :param isCartesian: coordinate type, either "lattice" - False or "cartesian" - True.
    :param threshold: threshold value for periodic boundary conditions

    :return: list of calculated bond lengths.
    """
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

    return delta_x, delta_y, delta_z


def plot_relative_coordinates(delta_x: List, delta_y: List, delta_z: List,
                              run_dir: str,
                              isCartesian: bool,
                              ymax: float, ymin: float,
                              show: bool,
                              dpi=300):
    """ Plot the bond lengths and save the plot to a file.

    :param delta_x: list of x-components.
    :param delta_y: list of y-components.
    :param delta_z: list of z-components.
    :param run_dir: directory where exciting runs.
    :param isCartesian: coordinate type, either "lattice" - False or "cartesian" - True.
    :param ymax:
    :param ymin:
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
        'lines.linewidth': 3.0,
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

    optimization_step = range(len(delta_x))

    plt.plot(delta_z, 'gd-', label=r'$\Delta$3')
    plt.plot(delta_y, 'bs-', label=r'$\Delta$2')
    plt.plot(delta_x, 'ro-', label=r'$\Delta$1')
    ax.set_xlabel('Optimization steps')
    ax.set_ylabel(f'Relative coordinate ({acoord})')
    ax.grid(True, linestyle='--')

    xmin = min(optimization_step)
    xmax = max(optimization_step)
    ymin_auto = min(min(delta_x), min(delta_y), min(delta_z))
    ymax_auto = max(max(delta_x), max(delta_y), max(delta_z))

    dxx = abs(xmax - xmin) / 18
    dyy = abs(ymax_auto - ymin_auto) / 15

    xmin = xmin - dxx
    xmax = xmax + dxx

    ymin = ymin if ymin is not None else ymin_auto - dyy
    ymax = ymax if ymax is not None else ymax_auto + dyy

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    ax.xaxis.set_major_locator(ptk.MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(ptk.ScalarFormatter(useOffset=True, useMathText=True))

    ax.legend()

    plt.savefig(join(run_dir, "PLOT.png"), format='png', dpi=dpi)

    if show:
        plt.show()

    plt.close()


def main() -> None:
    parser = ArgumentParser(
        description="Python visualization tool for following relative atomic coordinates of atoms during the relaxation process.")

    parser.add_argument("--run-directory", "-r",
                        default=os.getcwd(),
                        nargs=1,
                        dest="run_directory",
                        help="root path for files that are created by this script")

    parser.add_argument("--ATOM1", "--a1",
                        type=int,
                        help="First atom number",
                        default=1,
                        dest="atom1", )

    parser.add_argument("--ATOM2", "--a2",
                        type=int,
                        help="Second atom number",
                        default=2,
                        dest="atom2", )

    parser.add_argument("-sh", "--show",
                        action="store_true",
                        help="Show plot in a window")

    parser.add_argument("--y_min",
                        type=float,
                        dest="ymin",
                        help="minimum value of y axis",
                        default=None)

    parser.add_argument("--y_max",
                        type=float,
                        dest="ymax",
                        help="maximum value of y axis",
                        default=None)

    args = parser.parse_args()
    run_dir = args.run_directory
    atom1 = args.atom1
    atom2 = args.atom2
    ymax = args.ymax
    ymin = args.ymin
    show = args.show

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

    # Calculate relative coordinates
    delta_x, delta_y, delta_z = calculate_relative_coordinates(position_data, atom1, atom2, lattice_matrix, isCartesian)

    # Plot bond lengths
    plot_relative_coordinates(delta_x, delta_y, delta_z, run_dir, isCartesian, ymin, ymax, show)


if __name__ == "__main__":
    main()
