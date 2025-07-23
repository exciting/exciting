"""Fit energy-vs-volume curves.

Located at `excitingscripts/plot/volumecurves.py`.

Call as:

```bash
python3 -m excitingscripts.plot.volumecurves -r dir1 dir2 dir3
```
Where <code>dir1</code>, <code>dir2</code>, <code>dir3</code> take the place of the names of the  directories where
exciting calculations have been performed. The script can be used for any number of directories."""

import os
import pathlib
from argparse import ArgumentParser
from typing import Union

import matplotlib.pyplot as plt
import pylab as pyl

from excitingtools import ExcitingInputXML

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
ax1.set_ylabel(r'Energy - E$\mathregular{_{min}}$ [Ha]', labelpad=19)
ax1.set_xlabel(r'Volume [Bohr$\mathregular{^{3}}$]', labelpad=13)
pyl.grid(True)

def determine_functional(input_file: Union[str, pathlib.Path]) -> str:
    """Determine the name of the XC functional in a given input file.

    :param input_file: Input file.
    :returns: String containing the name of the XC functional.
    """

    parsed_input = ExcitingInputXML.from_xml(input_file)

    if hasattr(parsed_input.groundstate, "libxc"):
        try:
                functional = parsed_input.groundstate.libxc.correlation + "+" + parsed_input.groundstate.libxc.exchange
        except AttributeError:
                functional = parsed_input.groundstate.libxc.xc
    else:
        try:
            functional = parsed_input.groundstate.xctype
        except AttributeError:
            functional = "GGA_PBE_SOL"

    return functional

def plot_volumecurves(root_directory: Union[str, pathlib.Path], color_index) -> None:
    """Plot energy-vs-curve values for a given running directory.

    :param root_directory: Root directory.
    :param color_index: Index needing for plotting curves with different colors for each calculation.
    """
    functional = determine_functional(os.path.join(os.getcwd(), f"{root_directory}/input.xml"))

    with open(os.path.join(os.getcwd(), f"{root_directory}/energy-vs-volume"), "r") as f:
        lines = f.readlines()
    energy_volume_data = []

    for i in range(len(lines)):
        energy = (float(lines[i].split()[1]))
        strain = (float(lines[i].split()[0]))
        energy_volume_data.append([strain, energy])
    energy_volume_data = sorted(energy_volume_data)

    energy = [energy_volume_data[i][1] for i in range(len(energy_volume_data))]
    emin = min(energy)
    for i in range(len(energy)):
        energy[i] = energy[i] - emin
    strain = [energy_volume_data[i][0] for i in range(len(energy_volume_data))]

    ax1.plot(strain, energy, color=colors[color_index], label=functional, marker='o', markersize=12, linewidth=3.0,
             zorder=3)

def main() -> None:
    parser = ArgumentParser(description="Plot energy-vs-volume curves values for different functionals.")


    parser.add_argument("--root-directories", "-r",
                        nargs='+',
                        dest="root_directories",
                        help="names of root directories")

    parser.add_argument("-sh", "--show",
                        action="store_true",
                        help="show plot")

    args = parser.parse_args()

    for i in range(len(args.root_directories)):
        plot_volumecurves(args.root_directories[i], i)

    for line in ax1.get_xticklines() + ax1.get_yticklines():
        line.set_markersize(8)
        line.set_markeredgewidth(3)

    leg = ax1.legend(loc="upper center", borderaxespad=0.5, numpoints=1, fontsize=22)
    leg.get_frame().set_linewidth(4.0)
    leg.get_frame().set_edgecolor("grey")
    leg.draw_frame(True)

    plt.ylim(ymin=-0.001)
    xmin, xmax = plt.xlim()
    plt.hlines(0.0, xmin, xmax, linewidth=3.0, linestyles="dashed")

    plt.savefig('PLOT.png', orientation='portrait', format='png', dpi=300)

    if args.show:
        plt.show()

if __name__ == "__main__":
    main()
