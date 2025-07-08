"""Produce plot of the spin texture.

Located at `excitingscripts/plot/spintext.py`.

Call as:

```bash
python3 -m excitingscripts.plot.spintext -b ib -c context 
```
Where <code>ib</code> defines the band index for the plot and <code>context</code> defines the context of the contour plot. Choises are <code>energy</code> and <code>spin_z</code>.
"""

import os
import pathlib
from argparse import ArgumentParser
from typing import Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from excitingtools import ExcitingInputXML, parse
from matplotlib import gridspec


def triple_product(a: np.ndarray, b: np.ndarray, c: np.ndarray):
    r""" Vector triple product, defined as \mathbf{a} \cdot (\mathbf{b} \wedge \mathbf{c}).

    :param a: Vector a
    :param b: Vector b
    :param c: Vector c
    :return Triple product
    """
    return np.dot(a, np.cross(b, c))


def plane_transformation(rec_lat_vec: np.ndarray, plot_vec: np.ndarray):
    """Take reciprocal lattice vectors and ONS of a plane in rec. lat. coordinates where the first two vectors span the
    plane and the third is normal to them and calculate a matrix that transforms points in the plane to the xy plane
    in cartesian coordinates.

    :param rec_lat_vec: Reciprocal lattice vectors.
    :param plot_vec: ONS of the plotting plane.
    :return transformation_matrix: Matrix that transforms k and spin vectors to the plot plane.
    """
    norm = np.linalg.norm
    # transform plot vec in cartesian coordinates
    plot_vec = (rec_lat_vec.dot(plot_vec)).transpose()
    # extend plot vec to an orthogonal system
    plot_vec = np.array([(plot_vec[1] - plot_vec[0]) / norm(plot_vec[1] - plot_vec[0]),
                         (plot_vec[2] - plot_vec[0]) / norm(plot_vec[2] - plot_vec[0]),
                         np.cross(plot_vec[1] - plot_vec[0], plot_vec[2] - plot_vec[0]) \
                         / norm(np.cross(plot_vec[1] - plot_vec[0], plot_vec[2] - plot_vec[0]))])

    transformation_matrix = np.transpose(np.linalg.inv(plot_vec))

    return transformation_matrix


def reciprocal_lattice_vectors(lat_vec: np.ndarray):
    r"""Get the reciprocal lattice vectors of real-space lattice vectors \{\mathbf{a}\}:

    \mathbf{b}_0 = 2 \pi \frac{\mathbf{a}_1 \wedge \mathbf{a}_2} {\mathbf{a}_0 \cdot (\mathbf{a}_1 \wedge \mathbf{a}_2)}
    \mathbf{b}_1 = 2 \pi \frac{\mathbf{a}_2 \wedge \mathbf{a}_3} {\mathbf{a}_0 \cdot (\mathbf{a}_1 \wedge \mathbf{a}_2)}
    \mathbf{b}_2 = 2 \pi \frac{\mathbf{a}_0 \wedge \mathbf{a}_1} {\mathbf{a}_0 \cdot (\mathbf{a}_1 \wedge \mathbf{a}_2)}

    :param lat_vec: Lattice vectors, stored column-wise
    :return: rec_lat_vec: Reciprocal lattice vectors, stored column-wise
    """
    volume = triple_product(lat_vec[:, 0], lat_vec[:, 1], lat_vec[:, 2])
    rec_lat_vec = np.empty(shape=(3, 3))
    rec_lat_vec[:, 0] = 2 * np.pi * np.cross(lat_vec[:, 1], lat_vec[:, 2]) / volume
    rec_lat_vec[:, 1] = 2 * np.pi * np.cross(lat_vec[:, 2], lat_vec[:, 0]) / volume
    rec_lat_vec[:, 2] = 2 * np.pi * np.cross(lat_vec[:, 0], lat_vec[:, 1]) / volume

    return rec_lat_vec

def plot_spintext(root_directory: Union[str, pathlib.Path], band: int, contour: str,
                  contour_threshold: Union[str, float]):
    """ Plot the spin texture for a given band.

    :param root_directory: Directory of the exciting calculation.
    :param band:      Number of the band for the plot.
    :param contour:   Variable that will be plotted as contour. Can be either energy or spin_z.
    :param contour_threshold: Threshold for the contour plit. Can be either max or float.
    If max, the threshold is the absolute maximum value of the contour.
    """
    file_path_input = os.path.join(root_directory, 'input.xml')
    file_path_spintext = os.path.join(root_directory, 'spintext.xml')

    parsed_input = ExcitingInputXML.from_xml(file_path_input)

    if hasattr(parsed_input.structure.crystal_properties, "scale"):
        scale = parsed_input.structure.crystal_properties.scale
    else:
        scale = 1.0

    lat_vec = parsed_input.structure.lattice * scale

    parallelogram_input = parsed_input.properties.spintext.plot2d.parallelogram # pylint: disable=no-member
    plot_vec = np.array([parallelogram_input.origin.coord, parallelogram_input.point[0].coord,
                         parallelogram_input.point[1].coord]).transpose()
    grid = parallelogram_input.grid

    spin_text_data = parse(file_path_spintext)

    if spin_text_data[f"{len(spin_text_data) - 1}"]["ist"] >= band >= spin_text_data["0"]["ist"]:
        band_index = band - spin_text_data["0"]["ist"]
    else:
        raise ValueError('Band must be in the range of the bands considered for the spin texture.')

    rec_lat_vec = reciprocal_lattice_vectors(lat_vec)
    trans_mat = plane_transformation(rec_lat_vec, plot_vec)

    spin_text = spin_text_data[f"{band_index}"]
    k_point = [trans_mat.dot(np.array(k)) for k in spin_text["k-point"]]
    spin = [trans_mat.dot(np.array(s)) for s in spin_text["spin"]]
    energy = np.array(spin_text["energy"])

    bohr_to_angstrom = 0.529177

    k_x = np.array([k[0] / bohr_to_angstrom for k in k_point]).reshape(grid)
    k_y = np.array([k[1] / bohr_to_angstrom for k in k_point]).reshape(grid)

    s_x = np.array([s[0] for s in spin]).reshape(grid)
    s_y = np.array([s[1] for s in spin]).reshape(grid)
    s_z = np.array([s[2] for s in spin])

    mpl.rcParams['grid.linewidth'] = 3
    mpl.rcParams['xtick.labelsize'] = 25
    mpl.rcParams['ytick.labelsize'] = 25
    mpl.rcParams['xtick.major.width'] = 2
    mpl.rcParams['ytick.major.width'] = 2
    mpl.rcParams['xtick.major.size'] = 5
    mpl.rcParams['ytick.major.size'] = 5
    mpl.rcParams['axes.edgecolor'] = 'black'
    mpl.rcParams['axes.labelsize'] = 33  # fontsize of the x any y labels
    mpl.rcParams['axes.labelcolor'] = 'black'
    mpl.rcParams['axes.linewidth'] = 2.0  # set the value globally
    plt.rcParams['xtick.major.pad'] = 10
    plt.rcParams['ytick.major.pad'] = 10
    plt.rcParams.update({'mathtext.default': 'regular'})

    fig = plt.figure(figsize=(12, 10), dpi=300)
    gs = gridspec.GridSpec(1, 1, figure=fig)
    ax = fig.add_subplot(gs[0])

    # plot colour map
    if contour == "spin_z":
        if contour_threshold == "max":
            thr = round(max(np.abs(s_z)), 2)
        else:
            thr = float(contour_threshold)

        if thr == 0.0:
            raise ValueError(
                'Threshhold of contour plot for s_z is zero. You can set it manually with the -cthr option.')

        ticks = np.linspace(-1, 1, 5) * thr
        s_z = s_z.reshape(grid)

        cp = ax.contourf(k_x, k_y, s_z, 100, cmap="bwr", vmin=-thr, vmax=thr)

        cbar = fig.colorbar(cp, orientation="vertical", ticks=ticks)
        cbar.set_label('$s_z$')

    elif contour == "energy":
        thr_min, thr_max = round(min(energy), 1), round(max(energy), 1)
        ticks = thr_min + np.linspace(0.0, 1.0, 8) * (thr_max - thr_min)
        energy = np.array(energy).reshape(grid)
        cp = ax.contourf(k_x, k_y, energy, 100, cmap="autumn", vmin=thr_min, vmax=thr_max)
        cbar = fig.colorbar(cp, orientation="vertical", ticks=ticks)
        cbar.set_label('E [eV]', labelpad=50, rotation=-90)

    # plot spin texture
    ax.quiver(k_x, k_y, s_x, s_y)
    ax.set_xlabel("$k_1$ [a. u.]")
    ax.locator_params(nbins=3, axis='x')
    ax.set_ylabel("$k_2$ [a. u.]")
    ax.locator_params(nbins=3, axis='y')

    plt.savefig("PLOT-spintext.png", bbox_inches='tight')


def main() -> None:
    parser = ArgumentParser(description="Plot the spin texture for a given band.")

    parser.add_argument("-r",
                        type=Union[str, pathlib.Path],
                        default=[os.getcwd()],
                        dest="root_directory",
                        nargs=1,
                        help="path to the directory where input.xml and spintext.xml are stored")

    parser.add_argument('-b',
                        type=int,
                        default=0,
                        nargs=1,
                        dest="band",
                        help="band index for the plot")

    parser.add_argument('-c',
                        type=str,
                        nargs=1,
                        choices=['energy', 'spin_z'],
                        default=['spin_z'],
                        dest="contour",
                        help="Defines the context of the contour plot. Choices are energy and spin_z")

    parser.add_argument('-cthr',
                        type=str,
                        nargs=1,
                        default=['max'],
                        dest="contour_threshold",
                        help="threshhold for the contour in case of spin_z")

    parser.add_argument("-sh", "--show",
                        action="store_true",
                        help="show plot")

    args = parser.parse_args()

    try:
        cthr = float(args.contour_threshold[0])
    except:
        cthr = args.contour_threshold[0]
        if cthr != 'max':
            raise TypeError('contour_threshold needs to be float or "max"')

    plot_spintext(args.root_directory[0], args.band[0], args.contour[0], cthr)

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
