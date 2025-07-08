import os
from argparse import ArgumentParser
from os.path import join
from typing import Tuple, List, Dict, Callable

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ptk
import numpy as np
from excitingscripts.checkfit.energy_vs_strain import parse_info_elastic_constants
from excitingscripts.utils.utils import sort_lists_by_first_list, extract_values_from_line, birch_murnaghan_eos, \
    residuals, birch_murnaghan_fit
from scipy.constants import physical_constants, e

if matplotlib.__version__.split(".")[0] == "2":
    matplotlib.style.use("classic")

symmetry_factors = {0: 1,
                    1: 1,
                    2: 2,
                    3: 4}
symmetry_labels = {0: "     ",
                   1: "(sc) ",
                   2: "(bcc)",
                   3: "(fcc)"}

factor = 2

electron_charge = e
bohr_radius = physical_constants["Bohr radius"][0]
rydberg2ev = 13.605698066
unitconv = electron_charge * rydberg2ev / (1e9 * bohr_radius ** 3) * factor

# Format strings for output
fmt = '%10.4f'
amt = '%10.4f'
emt = '%9.5f'
bmt = '%8.3f'
pmt = '%16.10f'
lmt = '%10.2f'
a2t = '%12.3f'
a3t = '%14.3f'
al0 = '%7.4f'

# File names
volume_file = "energy-vs-volume"
strain_file = "energy-vs-strain"
planar_file = "planar"
source_file = "source.xml"
info_file = "INFO-elastic-constants"
output_file = "birch-murnaghan"

# Flags for existence of files
isStrain = os.path.exists(strain_file)
isPlanar = os.path.exists(planar_file)


def parse_energy_file(data_file: str) -> Tuple[List[float], List[float]]:
    """Reads the "energy-vs-volume" or "energy-vs-strain" file

    :param data_file: path to the file
    :return: strain/volume and energy values
    """
    if not os.path.exists(data_file):
        raise FileNotFoundError(data_file)

    X = []
    energy = []

    with open(data_file, "r") as input_energy:
        for line in input_energy:
            if len(line) != 0:
                x, ene = line.split()
                X.append(float(x))
                energy.append(float(ene))

    return X, energy

def strain_to_volume(strain, dim, v0):
    """

    :param strain: Strain Value
    :param dim: Dimension of the system
    :param v0: Volume at zero strain
    :return: Volume value
    """
    estrain = np.sqrt(1.0 + 2.0 * strain) - 1.0
    return (1 + estrain) ** dim * v0


def volume_to_strain(volume, dim, v0):
    """

    :param volume: Volume value
    :param dim: Dimension of the system
    :param v0: Volume at zero strain
    :return: Strain Value
    """
    estrain = (volume / v0) ** (1.0 / dim) - 1.0
    return ((estrain + 1.0) ** 2 - 1.0) / 2.0


def fit_energy_vs_volume(volumes: List[float], energies: List[float]) -> Tuple[
    Callable, List[float], List[float], List[float], Dict, float]:
    """
    Fit the energy-vs-volume or energy-vs-strain data using the Birch-Murnaghan Equation of State.

    :param volumes: List of volumes.
    :param energies: List of energies.
    :return:Birch fit, bulk modulus, bp, minima, lattice constant, chi value
    """
    if len(volumes) < 4:
        raise ValueError(f"\n ERROR: Too few volumes ({len(volumes)})!\n")

    p = birch_murnaghan_fit(volumes, energies)[0]
    v0, e0, b0, bp_ = p
    b0 = b0

    lattice_const = {i: [] for i in [1, 2, 3]}

    for isym in [1, 2, 3]:
        lattice_const[isym].append(np.cbrt(symmetry_factors[isym] * v0))

    chi = 0
    ebm = []
    fchi = 0
    for i in range(len(volumes)):
        chi = chi + residuals(p, energies[i], volumes[i]) ** 2
        fchi = fchi + (energies[i] - birch_murnaghan_eos(volumes[i], [v0, e0, b0 , bp_])) ** 2
        ebm.append(birch_murnaghan_eos(volumes[i], [v0, e0, b0, bp_]))

    chi = np.sqrt(chi) / len(volumes)

    def curv(xvol):
        fitted_energies = []
        for x in xvol:
            fitted_energies.append(birch_murnaghan_eos(x, p))
        return fitted_energies

    bulk_modulus = [b0 * unitconv]
    bp = [bp_]
    dmin = [v0]

    return curv, bulk_modulus, bp, dmin, lattice_const, chi


def print_info_to_stdout(dmin: List[float], lattice_const: Dict,
                         bulk_modulus: List[float], bp: List[float], chi: float, dim=3, v_eq: float = None) -> None:
    """Print information onto the screen

    :param dmin: minima of the fit.
    :param lattice_const: lattice constant for the given lattice symmetry code
    :param bulk_modulus: bulk modulus
    :param bp: derivative of bulk modulus with respect to pressure
    :param chi: chi-squared value indicating the goodness of fit
    :param dim:
    :param v_eq:
    """
    if len(dmin) > 1:
        print("##############################################\n")
        print("WARNING: Multiple minima are found!\n")
        print("##############################################\n")

    if len(dmin) == 0:
        print("##############################################\n")
        print("WARNING: No minimum in the given xrange!\n")
        print("##############################################\n")

    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    if isStrain:

        print("        A2            A3           lagrangian strain at minimum      log(chi)")

        bfactor = dim ** 2

        if isPlanar:
            with open(planar_file) as f:
                line = f.readline()

            alat, covera = extract_values_from_line(line)

        for i, v0 in enumerate(dmin):
            a2 = bulk_modulus[i] * bfactor
            a3 = -a2 * (dim * (bp[i] - 2.0) + 6.0)
            mls = volume_to_strain(v0, dim, v_eq)
            if isPlanar:
                print(
                    f"{a2t % a2} {a3t % a3}      {emt % mls}  ( alat0 =  {al0 % ((1.0 + mls) * alat)} ) {lmt % (np.log10(chi))}")
            else:
                print(f"{a2t % a2} {a3t % a3}                {emt % mls}            {lmt % (np.log10(chi))}")

    else:

        print("     V0        B0         Bp        a-sc       a-bcc      a-fcc     log(chi)")
        for i, v0 in enumerate(dmin):
            print(fmt % v0, bmt % (bulk_modulus[i]), bmt % (bp[i]), amt % (lattice_const[1][i]),
                  amt % (lattice_const[2][i]), amt % (lattice_const[3][i]), lmt % (np.log10(chi)))

    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    with open(output_file, 'w') as output:
        for i, v0 in enumerate(dmin):
            output.write(pmt % v0 + ' ' + pmt % bulk_modulus[i] + ' ' + pmt % bp[i] + ' ')
            if isStrain:
                mls = volume_to_strain(v0, dim, v_eq)
                output.write(pmt % mls + ' ' + pmt % mls + ' ' + pmt % mls + ' ')
            else:
                output.write(
                    pmt % lattice_const[1][i] + ' ' + pmt % lattice_const[2][i] + ' ' + pmt % lattice_const[3][
                        i] + ' ')
            output.write(pmt % np.log10(chi) + '\n')


def plot_energy(X: List[float], energies: List[float], curv: Callable, dmin: np.ndarray,
                output_dir: str, dim: int, v_eq: float):
    """
    Plot the energy-vs-volume or energy-vs-strain data and save the plot.

    :param X : List of volumes.
    :param energies: List of energies
    :param curv : Fitted curve.
    :param dmin : Minima volume of the fit.
    :param output_dir : Directory to save the plot.
    :param dim: Dimension of the system
    :param v_eq: Volume at zero strain
    """
    plt.rcParams.update({
        'figure.figsize': (10, 7.5),
        'ytick.minor.size': 6,
        'xtick.major.pad': 8,
        'ytick.major.pad': 4,
        'patch.linewidth': 2.,
        'axes.linewidth': 2.,
        'lines.linewidth': 1.8,
        'lines.markersize': 8.0,
        'axes.formatter.limits': (-4, 6),
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'axes.labelsize': 20,
        'axes.labelcolor': 'black',
        'axes.axisbelow': 'True',
        'legend.fontsize': 22,
        'axes.titlesize': 24,
        'axes.titlepad': 20
    })
    # Create figure and axis
    fig, ax = plt.subplots()

    yfmt = ptk.ScalarFormatter(useOffset=True, useMathText=True)

    xvalues = np.linspace(min(X), max(X), 100)
    energy_values = curv(xvalues)
    energy_min = curv(dmin)

    if isStrain:
        ax.set_xlabel(u'Lagrangian strain')

        # Convert Volume values to Strain values
        for i in range(len(dmin)):
            dmin[i] = volume_to_strain(dmin[i], dim, v_eq)
        for i in range(len(X)):
            X[i] = volume_to_strain(X[i], dim, v_eq)
        for i in range(len(xvalues)):
            xvalues[i] = volume_to_strain(xvalues[i], dim, v_eq)

    else:
        ax.set_xlabel(u'Volume [Bohr\u00B3]')

    ax.set_ylabel(r'Energy [Ha]')

    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(2)

    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.grid(True)

    plt.plot(xvalues, energy_values, 'b-', label=f'birch-murnaghan fit')
    plt.plot(X, energies, 'go', label='calculated')
    plt.plot(dmin, energy_min, 'ro')
    plt.legend(loc=9, borderaxespad=.8, numpoints=1)

    ymax = max(max(energy_values), max(energies))
    ymin = min(min(energy_values), min(energies))
    dxx = abs(max(xvalues) - min(xvalues)) / 18
    dyy = abs(ymax - ymin) / 18
    ax.yaxis.set_major_formatter(yfmt)
    ax.set_xlim(min(xvalues) - dxx, max(xvalues) + dxx)
    ax.set_ylim(ymin - dyy, ymax + dyy)

    ax.xaxis.set_major_locator(ptk.MaxNLocator(7))

    plt.tight_layout()
    plt.savefig(join(output_dir, 'PLOT.png'), orientation='portrait', format='png', dpi=300)

    plt.close()


def main() -> None:
    parser = ArgumentParser(
        description="Python script for fitting energy-vs-volume curves using BM-EoS")

    parser.add_argument("--root-directory", '-r',
                        type=str,
                        default=[os.getcwd()],
                        nargs=1,
                        dest="directory",
                        help="Directory containing the necessary files")

    args = parser.parse_args()
    directory = args.directory[0]

    volume_file = join(directory, "energy-vs-volume")
    strain_file = join(directory, "energy-vs-strain")
    planar_file = join(directory, "planar")
    isStrain = os.path.exists(strain_file)
    isPlanar = os.path.exists(planar_file)

    if not isStrain:
        data_file = volume_file
        print("\n The input file is \"energy-vs-volume\".")
    else:
        data_file = strain_file
        print("\n The input file is \"energy-vs-strain\".")

        info = parse_info_elastic_constants(directory)

        deformation_code = int(info["Deformation code"])
        v_eq = info["Volume of equilibrium unit cell"]

        if deformation_code in [1, 2, 3]:
            dim = 1
        elif deformation_code == 8:
            dim = 2
        elif deformation_code == 0:
            dim = 3
        else:
            print("\n ERROR: deformation type " + str(deformation_code) + " not allowed!\n")
            return
        if isPlanar:
            print(" Modified version for strained planar systems!\n")

    # Read data from energy file
    X, energies = parse_energy_file(data_file)

    # Converting strain to volume
    volumes = []
    for i, x in enumerate(X):
        if not isStrain:
            volumes.append(x)
        else:
            vol = strain_to_volume(x, dim, v_eq)
            volumes.append(vol)

    if isPlanar:
        with open(planar_file) as f:
            line = f.readline()

        alat, covera = extract_values_from_line(line)
        global unitconv
        unitconv = unitconv * alat * covera * 5.2917721092e-2
    else:
        dim = 3
        v_eq = None

    volumes, energies = sort_lists_by_first_list(volumes, energies)
    curv, bulk_modulus, bp, dmin, lattice_const, chi = fit_energy_vs_volume(volumes, energies)
    print_info_to_stdout(dmin, lattice_const, bulk_modulus, bp, chi, dim, v_eq)
    plot_energy(volumes, energies, curv, dmin, directory, dim, v_eq)


if __name__ == "__main__":
    main()
