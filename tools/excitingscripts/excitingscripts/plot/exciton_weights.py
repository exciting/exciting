"""Visualize energy-vs-strain curves.

Located at `excitingscripts/plot/exciton_weights.py`. 

Call as:

```bash
python3 -m excitingscripts.plot.exciton_weights structure_name file_name energy_min energy_max exciton_weights_size
```
Where <code>structure_name</code> is the name of the structure, <code>file_name</code> is the name of the file containing data needed for exciton visualization, <code>energy_min</code> and <code>energy_max</code> are the minimum and maximum energy values for setting plot axis limits, and <code>exciton_weights_size</code> is the size of excitonic weights.
"""

import pathlib
from argparse import ArgumentParser
from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
from excitingtools.exciting_obj_parsers.ks_band_structure import parse_band_structure
from scipy.constants import physical_constants


def plot_exciton_weights(structure_name: str, exciton_file : Union[str, pathlib.Path], energy_min: float,
                         energy_max: float, exciton_weights_size: float) -> None:
    """Plot excitonic weights along a band structure path.

    Assumes presence of "bandstructure.dat" file in current running directory for plotting band structure.

    :param structure_name: Name of structure.
    :param exciton_file: File containing data needed for exciton visualization.
    :param energy_min: Minimum energy value for setting plot axis limit.
    :param energy_max: Maximum energy value for setting plot axis limit.
    :param exciton_weights_size: Size of excitonic weights.
    """
    band_structure_file = Path("bandstructure.dat")
    if not band_structure_file.exists():
        raise FileNotFoundError("bandstructure.dat file not found")

    if not Path(exciton_file).exists():
        raise FileNotFoundError(f"{exciton_file} file not found")

    ha_to_ev = physical_constants["hartree-electron volt relationship"][0]

    with open(band_structure_file) as f:
        header = f.readline()

    number_kpts = int(header.split()[3])
    number_bands = int(header.split()[2])

    band_data = parse_band_structure(band_structure_file)
    band_data.bands = band_data.bands * ha_to_ev
    band_data.e_fermi = 0.0

    vertices, labels = band_data.band_path()

    exciton_weights_data = np.genfromtxt(exciton_file)
    exciton_weights = np.reshape(exciton_weights_data[:, 2], (number_kpts, number_bands), order='F')

    # Plot settings
    figcolor = "white"
    fig = plt.figure(figsize=(10, 10))
    fig.patch.set_edgecolor(figcolor)
    fig.patch.set_facecolor(figcolor)

    plt.rcParams['axes.linewidth'] = 3.0
    plt.rcParams['grid.linewidth'] = 1.5
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30
    plt.rcParams['axes.edgecolor'] = 'black'
    plt.rcParams['axes.labelsize'] = 30
    plt.rcParams['axes.labelcolor'] = 'black'
    plt.rcParams['axes.axisbelow'] = 'True'
    plt.rcParams['legend.fontsize'] = 25
    plt.rcParams['xtick.major.pad'] = 10
    plt.rcParams['ytick.major.pad'] = 10
    plt.rcParams['axes.titlesize'] = 30
    plt.rcParams['axes.titlepad'] = 20

    ax1 = fig.add_axes([0.17, 0.1, 0.75, 0.8])
    ax2 = ax1.twinx()

    ax1.xaxis.grid(True, which="major", color="k", linestyle="-", linewidth=2)
    ax1.xaxis.set_label_position("bottom")
    ax1.set_xlim([vertices[0], vertices[-1]])
    ax1.set_xticks(ticks=vertices, labels=labels)
    ax1.set_ylim([energy_min, energy_max])
    ax1.set_ylabel('Energy [eV]')

    scale_factor = 200
    scale = exciton_weights_size * scale_factor

    # Band structure and excitonic weights plot
    ax1.plot(band_data.flattened_k_points, band_data.bands, color="b", lw=3.0, zorder=10)
    for i in range(number_bands):
        ax1.scatter(band_data.flattened_k_points, band_data.bands[:, i], s=(scale * exciton_weights[:, i]) ** 2, lw=3.0,
                    edgecolor='r', facecolor='none', zorder=11)

    ax1.grid(True)
    plt.title(f"{structure_name} excitonic weights")

    for line in ax1.get_xticklines() + ax1.get_yticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(2)

    # Fermi energy
    ax1.plot([vertices[0], vertices[-1]], [band_data.e_fermi, band_data.e_fermi], 'k', lw=3.0, ls='-')
    ax2.set_ylim(ax1.get_ylim())
    ax2.set_yticks([band_data.e_fermi], labels=[r'$\mathregular{E_F}$'])

def main() -> None:
    parser = ArgumentParser(description="Plot excitonic weights along a band structure path.")


    parser.add_argument("structure_name",
                        type=str,
                        nargs=1,
                        help="structure name")

    parser.add_argument("exciton_file",
                        nargs=1,
                        help="file containing data needed for exciton visualization")

    parser.add_argument("energy_min",
                        type=float,
                        nargs=1,
                        help="energy minimum")

    parser.add_argument("energy_max",
                        type=float,
                        nargs=1,
                        help="energy maximum")

    parser.add_argument("exciton_weights_size",
                        type=float,
                        nargs=1,
                        help="size of excitonic weights")

    parser.add_argument("-sh", "--show",
                        action="store_true",
                        help="show plot")

    args = parser.parse_args()

    plot_exciton_weights(args.structure_name[0], args.exciton_file[0], args.energy_min[0], args.energy_max[0],
                         args.exciton_weights_size[0])

    plt.savefig(f'{args.exciton_file[0]}.png', orientation='portrait', format='png', dpi=300)

    if args.show:
        plt.show()

if __name__ == "__main__":
    main()
