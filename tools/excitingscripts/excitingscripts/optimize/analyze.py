"""Python script for fitting the energy-vs-volume and energy-vs-strain curves."""

import warnings
from argparse import ArgumentParser
from pathlib import Path
from typing import List, Callable

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ptk
import numpy as np
from excitingscripts.optimize.setup import (
    get_deformation_matrix,
    optimization_map,
)
from excitingscripts.utils.utils import (
    birch_murnaghan_eos,
    pressure_birch_murnaghan_eos,
    initial_guess,
)
from excitingscripts.utils.utils import sort_lists_by_first_list
from excitingtools import ExcitingInputXML
from excitingtools import parse
from scipy.constants import physical_constants
from scipy.optimize import leastsq

if matplotlib.__version__.split(".")[0] == "2":
    matplotlib.style.use("classic")

bohr_radius = physical_constants["Bohr radius"][0]
joule2hartree = physical_constants["hartree-joule relationship"][0]
unitconv = joule2hartree / (bohr_radius ** 3) * 10 ** -9


def get_energy_data(run_dir: Path, dir_name: str, vol_flag: bool) -> None:
    """Extracts and saves the energy-volume data from a series of exciting calculations.

    :param run_dir: Path to the directory containing the calculations.
    :param dir_name: Name of the directory indicating the type of optimization.
    :param vol_flag: Flag to indicate if the calculation is a volume optimization.
    """
    displ_points = sum(
        1
        for directory in run_dir.iterdir()
        if directory.is_dir() and directory.name.startswith(f"{dir_name}_")
    )

    total_energy = []
    x_values = []

    if vol_flag:
        outfile = run_dir / "energy-vs-volume"
    else:
        outfile = run_dir / "energy-vs-strain"

    for i in range(displ_points):
        root_directory = run_dir / f"{dir_name}_{i + 1}"
        infofile = root_directory / "INFO.OUT"
        parsed_info = parse(infofile.as_posix())
        max_scf = max([int(i) for i in parsed_info["scl"].keys()])
        converged_results = parsed_info["scl"][str(max_scf)]
        total_energy.append(converged_results["Total energy"])

        if vol_flag:
            x_values.append(parsed_info["initialization"]["Unit cell volume"])

    if not vol_flag:
        strain_infofile = run_dir / f"{dir_name}-Parameters"
        lines = strain_infofile.read_text().splitlines()
        for line in lines:
            if "Physical strain" in line:
                x_values.append(float(line.split()[-1]))

    energy_data = np.vstack((x_values, total_energy)).T

    with open(outfile, "w") as f:
        np.savetxt(f, energy_data)


def murnaghan_eos(v: float, p: List[float]) -> float:
    """Calculate energy using the Murnaghan equation of state.

    :param v: Volume.
    :param p: Parameters for the equation [v0, e0, b0, bp].
    :return: Calculated energy.
    """
    eq_vol = p[0]  # Equilibrium volume
    min_energy = p[1]  # Minimum energy
    bulk_modulus = p[2]  # Bulk modulus at equilibrium volume
    bulk_modulus_pressure_deriv = p[3]  # Pressure derivative of bulk modulus

    energy = (
            min_energy
            + bulk_modulus
            * v
            / bulk_modulus_pressure_deriv
            * (
                    1
                    / (bulk_modulus_pressure_deriv - 1)
                    * (eq_vol / v) ** bulk_modulus_pressure_deriv
                    + 1
            )
            - bulk_modulus * eq_vol / (bulk_modulus_pressure_deriv - 1)
    )
    return energy


def pressure_murnaghan_eos(v: float, p: List[float]) -> float:
    """Calculate pressure using the Murnaghan equation of state.

    :param v: Volume.
    :param p: Parameters for the equation [v0, e0, b0, bp].
    :return: Calculated pressure.
    """
    eq_vol = p[0]  # Equilibrium volume
    bulk_modulus = p[2]  # Bulk modulus at equilibrium volume
    bulk_modulus_pressure_deriv = p[3]  # Pressure derivative of bulk modulus

    pressure = (
            3
            / 2
            * bulk_modulus
            * ((eq_vol / v) ** (7 / 3) - (eq_vol / v) ** (5 / 3))
            * (
                    1
                    + 3 / 4 * (bulk_modulus_pressure_deriv - 4) * ((eq_vol / v) ** (2 / 3) - 1)
            )
    )
    return pressure


def residuals(p: List[float], e: float, v: float, eos: Callable) -> float:
    """Calculate residuals for the least squares fit.

    :param p: Parameters for the Murnaghan EOS.
    :param e: Energies.
    :param v: Volumes.
    :param eos: Equation of state.
    :return: Residuals.
    """
    # Ensure e and v are numpy arrays for element-wise operations
    e = np.asarray(e)
    v = np.asarray(v)

    # Calculate residuals for each volume
    residual = np.array(
        [e_single - eos(v_single, p) for e_single, v_single in zip(e, v)]
    )

    return residual


def fit_energy_vs_volume(
        volumes: List[float], energies: List[float], out_file: Path, eos: str
) -> List:
    """Fit the energy data to an equation of state and return the optimized structure.

    :param volumes: List of volumes.
    :param energies: List of energies.
    :param out_file: Path to the output file.
    :param eos: Equation of state.
    :return: Fitted parameters.
    """
    if eos == "m":
        eos = murnaghan_eos
        pressure_eos = pressure_murnaghan_eos
        eos_string = " === Murnaghan eos ==============================="

    else:
        eos = birch_murnaghan_eos
        pressure_eos = pressure_birch_murnaghan_eos
        eos_string = " === Birch-Murnaghan eos ========================="

    # Initial guess
    p0 = initial_guess(volumes, energies)

    # Perform the least squares fit
    p, cov_x, infodict, mesg, ier = leastsq(
        residuals, p0, args=(energies, volumes, eos), full_output=True
    )

    eq_vol, min_energy, bulk_modulus, bulk_modulus_pressure_deriv = p

    # Calculate the residue
    residue = sum(infodict["fvec"] ** 2)

    # Calculate fit accuracy
    fit_accuracy = np.log10(np.sqrt(residue))

    bulk_modulus = bulk_modulus * unitconv  # Convert B0 to GPa units

    output_str = (
        " Fit accuracy:\n"
        f"     Log(Final residue in [Ha]): {round(fit_accuracy, 2)}\n\n"
        " Final parameters:\n"
        f"     E_min =  {round(min_energy, 7)} [Ha]\n"
        f"     V_min =  {round(eq_vol, 3)} [Bohr^3]\n"
        f"     B_0   =  {round(bulk_modulus, 3)} [GPa]\n"
        f"     B'    =  {round(bulk_modulus_pressure_deriv, 3)}\n"
    )

    print(output_str)

    with open(out_file, "w") as f:
        f.write(eos_string + "\n")
        f.write(output_str)
        f.write(" =================================================\n")
        f.write("\n Volume     E_dft-E_eos     Pressure [GPa]\n")
        for i, v in enumerate(volumes):
            e_diff = energies[i] - eos(v, p)
            pressure = pressure_eos(v, p)
            f.write(f"{v:10.4f} {e_diff:+13.8f} {pressure:10.3f}\n")

    def curv(xvol):
        fitted_energies = []

        # Check if xvol is iterable
        if not hasattr(xvol, "__iter__"):
            xvol = [xvol]

        for x in xvol:
            fitted_energy = eos(x, p)
            fitted_energies.append(fitted_energy)

        if len(fitted_energies) == 1:
            fitted_energies = fitted_energies[0]
        return fitted_energies

    return curv, bulk_modulus, bulk_modulus_pressure_deriv, eq_vol, fit_accuracy


def fit_energy_vs_strain(
        strain: List[float], energies: List[float], out_file: Path, order_of_fit=4
) -> List[float]:
    """Fit the energy data to a polynomial and return the optimized structure.

    :param strain: List of strains.
    :param energies: List of energies.
    :param out_file: Path to the output file.
    :param order_of_fit: Order of the polynomial fit.
    :return: Fitted parameters
    """
    # Check if the order of the fit is less than the number of data points
    if order_of_fit >= len(strain):
        raise ValueError(
            f"The order of polynomial({order_of_fit}) must be less than the number of data points({len(strain)})."
        )

    # Fit the data
    coefficients = np.polyfit(strain, energies, order_of_fit)
    curv = np.poly1d(coefficients)

    s_fit = np.linspace(strain[-1] * -1.2, strain[-1] * 1.2, 1000)
    e_fit = curv(s_fit)
    s0 = s_fit[e_fit.argmin()]

    return curv, None, None, s0, None


def fit_energy_data(run_dir: Path, vol_flag, dir_name: str, fit_type=None) -> None:
    """Fit the energy data and write the optimized structure to a new .xml file.

    :param run_dir: Path to the directory containing the calculations.
    :param vol_flag: Flag to indicate if the calculation is a volume optimization.
    :param dir_name: Name of the directory indicating the type of optimization.
    :param fit_type: Type of fit.
    """
    # Check for the equation of state
    if vol_flag:
        if fit_type is None or fit_type.lower() not in ("m", "b"):
            raise ValueError(
                "The equation of state must be provided for volume optimizations. Choose from Murnaghan (M/m) or "
                "Birch-Murnaghan (B/b) equation of state"
            )
        fit_type = fit_type.lower()
    else:
        if fit_type is None:
            fit_type = 4
        elif fit_type is not int:
            warnings.warn("The order of the polynomial fit must be an integer.")
            fit_type = 4

    # Load the energy data
    if vol_flag:
        data_file = run_dir / "energy-vs-volume"
        fit = fit_energy_vs_volume
    else:
        data_file = run_dir / "energy-vs-strain"
        fit = fit_energy_vs_strain

    X = []
    energy = []

    lines = data_file.read_text().splitlines()
    for line in lines:
        if len(line) != 0:
            x, ene = line.split()
            X.append(float(x))
            energy.append(float(ene))

    # Sort the data
    X, energy = sort_lists_by_first_list(X, energy)

    print("\n =====================================================================")

    out_name = dir_name.lower()

    if vol_flag:
        if fit_type.lower() == "m":
            eq_of_state = "M"
        else:
            eq_of_state = "BM"
        out_name = f"{eq_of_state}_eos"

    out_file = run_dir / f"{out_name}.out"

    # Fit the data. The parameters other than will be none for strain fit
    curv, b0, bp, xmin, fit_accuracy = fit(X, energy, out_file, fit_type)

    # Write the optimized .xml file
    # Parse the source.xml file
    source_file = run_dir / "source.xml"

    # Check if the source file exists
    if not source_file.exists():
        raise FileNotFoundError(f"Source file {source_file} not found.")

    # Extract base vectors from the input file
    parsed_source = ExcitingInputXML.from_xml(source_file)

    # Extract crystal properties
    if hasattr(parsed_source.structure.crystal_properties, "scale"):
        ref_scale = parsed_source.structure.crystal_properties.scale
    else:
        ref_scale = 1.0

    if hasattr(parsed_source.structure.crystal_properties, "stretch"):
        stretch = parsed_source.structure.crystal_properties.stretch
    else:
        stretch = [1.0, 1.0, 1.0]

    base_vectors = np.array(parsed_source.structure.lattice)

    # Calculate initial volume
    v_in = abs(np.linalg.det(base_vectors) * ref_scale ** 3 * np.prod(stretch))

    if vol_flag:
        # Calculate the strain
        strain = np.cbrt(xmin / v_in) - 1.0
    else:
        strain = xmin

    transformation_matrix = get_deformation_matrix(strain)[dir_name.upper()]

    # Update the base vectors
    new_base_vectors = np.dot(base_vectors, transformation_matrix)

    # Update base vectors in the XML
    parsed_source.structure.lattice = new_base_vectors.tolist()

    # Write the distorted structure to the XML file
    if vol_flag: out_name = eq_of_state
    parsed_source.write(run_dir / f"{out_name}-optimized.xml")
    print(
        f' Optimized lattice parameter saved into the file: "{out_name}-optimized.xml"'
    )
    print(" =====================================================================\n")

    # Plot the data
    plot_energy(energy, X, run_dir, curv, [b0, bp, xmin], vol_flag, out_name, fit_type)


def plot_energy(
        energies: List[float],
        x_data: List[float],
        output_dir: Path,
        curv: Callable,
        p: List[float],
        vol_flag: bool,
        out_name: str,
        fit_type,
) -> None:
    """Plot the fitted energy data and save the plot

    :param energies: List of energies.
    :param x_data: List of volumes/strain.
    :param output_dir: Directory to save the plot.
    :param curv: Fitted curve.
    :param p: Fitted parameters.
    :param vol_flag: Flag to indicate if the calculation is a volume optimization.
    :param out_name:
    :param fit_type: Type of fit.
    """
    plt.rcParams.update(
        {
            "figure.figsize": (10, 7.5),
            "ytick.minor.size": 6,
            "xtick.major.pad": 8,
            "ytick.major.pad": 4,
            "patch.linewidth": 2.0,
            "axes.linewidth": 2.0,
            "lines.linewidth": 3.5,
            "lines.markersize": 12.0,
            "axes.formatter.limits": (-4, 6),
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "axes.labelsize": 20,
            "axes.labelcolor": "black",
            "axes.axisbelow": "True",
            "legend.fontsize": 22,
            "axes.titlesize": 24,
            "axes.titlepad": 20,
            "font.size": 22,
        }
    )
    # Create figure and axis
    fig, ax = plt.subplots()

    yfmt = ptk.ScalarFormatter(useOffset=True, useMathText=True)

    e_min = curv(p[-1])

    if vol_flag:
        xlabel = "Volume [Bohr\u00b3]"
        fit_label = "Murnaghan eos" if fit_type == "m" else "Birch-Murnaghan eos"
        plot_text = (
            "E$_{{min}}$ = {:.7f} Ha\n"
            "V$_{{min}}$ = {:.3f} Bohr$^3$\n"
            "B$_0$ = {:.3f} GPa\n"
            r"B$^\prime$ = {:.5f}"
        ).format(e_min, p[2], p[0], p[1])
    else:
        xlabel = "Physical strain $\\epsilon$"
        fit_label = f"Order {fit_type} polynomial fit"
        plot_text = ("E$_{{min}}$ = {:.7f} Ha\n" "$\\epsilon_{{min}} = {:.5f}$").format(
            e_min, p[2]
        )

    ylabel = "Energy [Ha]"

    # Set labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Set the grid
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(2)

    plt.grid(True, linestyle="--")

    # Generate the data for the plot
    x_mn = min(x_data)
    x_mx = max(x_data)
    dx = (x_mx - x_mn) * 0.1

    xx = np.linspace(x_mn - dx, x_mx + dx, 1000)
    yy = curv(xx)

    x0 = x_data
    y0 = energies

    # Plot the fitted curve
    ax.plot(xx, yy, color="red", linewidth=2, label=fit_label)

    # Plot the DFT calculated points
    ax.plot(
        x0,
        y0,
        "o",
        color="green",
        markersize=8,
        markeredgecolor="black",
        markeredgewidth=1,
        label="DFT Calc.",
    )

    # Set legend
    ax.legend(numpoints=1, loc=9)

    # Set the text
    # Keep it below the legend box
    plt.text(
        0.3,
        0.75,
        plot_text,
        transform=ax.transAxes,
        verticalalignment="top",
    )

    # Set the plot limits
    max_y = max(max(yy), max(y0))
    min_y = min(min(yy), min(y0))

    max_x = max(max(xx), max(x0))
    min_x = min(min(xx), min(x0))

    dyy = (max_y - min_y) / 15
    ax.set_ylim(min_y - dyy, max_y + dyy)
    dxx = (max_x - min_x) / 18

    ax.set_xlim(min_x - dxx, max_x + dxx)

    ax.yaxis.set_major_formatter(yfmt)
    ax.xaxis.set_major_locator(ptk.MaxNLocator(7))

    out_name = out_name + ".png"

    plt.tight_layout()
    plt.savefig(output_dir / out_name, orientation="portrait", format="png", dpi=300)

    plt.close()


def main() -> None:
    parser = ArgumentParser(
        description="Python script for fitting the energy-vs-volume and energy-vs-strain curves."
    )

    parser.add_argument(
        "-d",
        "--directory",
        type=str,
        default=Path.cwd(),
        help="Directory where exciting runs",
    )

    parser.add_argument(
        "--eos",
        "--equation-of-state",
        type=str,
        default=None,
        help="""The equation of state to be used for fitting the data. Required for volume optimizations.
              Choose from Murnaghan (M/m) or Birch-Murnaghan (B/b) equation of state.
              If not provided, the script will default to a 4th order polynomial fit for strain optimizations.""",
    )

    args = parser.parse_args()
    run_dir = Path(args.directory)
    eos = args.eos

    fit_type = None  # integer: order of polynomial for energy-vs-strain fit or string: eos for energy-vs-volume fit

    # Check the type of optimization based on directory name
    vol_flag = False  # Is the calculation a volume optimization or other kind?

    dir_name = run_dir.stem.lower()
    optimization_keys = optimization_map.keys()

    if "vol" in dir_name.lower():
        vol_flag = True
        dir_name = "vol"
        try:
            fit_type = eos.lower()
        except AttributeError:
            warnings.warn(
                "The equation of state must be provided for volume optimizations. Defaulting to Murnaghan."
            )
            fit_type = "m"
    else:
        for key in optimization_keys:
            if key.lower() in dir_name and key != "VOL":
                dir_name = (
                    key.lower()
                )  # Set the directory name to the optimization type
                break
        try:
            fit_type = int(eos)
        except (ValueError, TypeError):
            pass

    # Extract the energy data
    get_energy_data(run_dir, dir_name, vol_flag)

    # Fit the data
    fit_energy_data(run_dir, vol_flag, dir_name, fit_type)


if __name__ == "__main__":
    main()
