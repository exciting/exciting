"""Visualize multiple energy-vs-distance curves.

Located at `excitingscripts/plot/compare_vdW.py`.

Call as:

```bash
python3 -m excitingscripts.plot.compare_vdW -f file_name -r dir1 dir2 dir3
```
"""

import os
import pathlib
from argparse import ArgumentParser
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pylab as pyl

# Plot settings
figcolor = 'white'
dpi = 300
fig = plt.figure(figsize=(14.5,10),dpi=dpi)
fig.patch.set_edgecolor(figcolor)
fig.patch.set_facecolor(figcolor)

plt.rcParams['axes.linewidth' ] = 4.0
plt.rcParams['grid.linewidth' ] = 1.5
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30
plt.rcParams['axes.edgecolor' ] = 'black'
plt.rcParams['axes.labelsize' ] = 45
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['axes.axisbelow' ] = 'True'
plt.rcParams['legend.fontsize'] = 30
plt.rcParams['xtick.major.pad'] = 20
plt.rcParams['ytick.major.pad'] = 10

colors=['b','r','g','y','k']

ax1 = fig.add_axes([0.2,0.18,0.75,0.76])

ax1.xaxis.set_label_position('bottom')
ax1.set_ylabel('Energy [Ha]', labelpad=19)
ax1.set_xlabel('Interlayer distance [Bohr]', labelpad=13)
pyl.grid(True)


def plot_compare_vdw(plot_file_path : Union[str, pathlib.Path], color_index) -> None:
    """Plot binding energy curve values for a given running directory.

    :param plot_file_path: Path to file containing data wanted for plot.
    :param color_index: Index needing for plotting curves with different colors for each calculation.
    """
    energy_strain_data = np.genfromtxt(plot_file_path)
    root_directory = os.path.dirname(plot_file_path)
    energy_strain_data = energy_strain_data[energy_strain_data[:, 0].argsort()]
    ax1.plot(energy_strain_data[:, 0], energy_strain_data[:, 1], color=colors[color_index], label=root_directory,
             marker='o', markersize=12, linewidth=3.0, zorder=3)

def main() -> None:
    parser = ArgumentParser(description="Plot binding energy curves values for different van-der-Waals Corrections .")

    parser.add_argument("--root-directories", "-r",
                        nargs='+',
                        dest="root_directories",
                        help="names of root directories")

    parser.add_argument("--plot-file", "-f",
                        nargs=1,
                        dest="plotfile",
                        help="name of the file containg plot data")

    parser.add_argument("-sh", "--show",
                        action="store_true",
                        help="show plot")

    args = parser.parse_args()

    for index, directory in enumerate(args.root_directories):
        plot_compare_vdw(os.path.join(directory, args.plotfile[0]), index)

    for line in ax1.get_xticklines() + ax1.get_yticklines():
        line.set_markersize(8)
        line.set_markeredgewidth(3)

    leg = ax1.legend(loc="upper center", borderaxespad=0.5, numpoints=1, fontsize=22)
    leg.get_frame().set_linewidth(4.0)
    leg.get_frame().set_edgecolor("grey")
    leg.draw_frame(True)

    plt.savefig('PLOT.png', orientation='portrait', format='png', dpi=300)

    if args.show:
        plt.show()

if __name__ == "__main__":
    main()
