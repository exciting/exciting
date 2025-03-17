from argparse import ArgumentParser
from typing import Tuple

import numpy as np
from scipy.constants import physical_constants

# Constants
inv_cm_2_ev = physical_constants["inverse meter-electron volt relationship"][0] * 10**2
k_2_ev = physical_constants["kelvin-electron volt relationship"][0]
ev_2_ha = physical_constants["electron volt-hartree relationship"][0]
ha_2_ev = physical_constants["hartree-electron volt relationship"][0]
ha_2_inv_cm = physical_constants["hartree-inverse meter relationship"][0] * 10**(-2)
eps_min = 1.e-15
eps_max = 50

def heat_capacity_distribution(omega: float) -> float:
    """Calculate distribution which is later multiplied by the phonon DOS to determine the heat capacity.

    :param omega: Frequency.
    :returns Distribution which is needed for determining the heat capacity.
    """
    omega_dist = 1
    if omega > eps_min:
        omega_dist = 0
        if omega < eps_max:
            omega_dist = omega**2 * np.exp(omega) / (np.exp(omega) - 1)**2

    return omega_dist

def normalized_frequency(temp: float, omega: float)-> float:
    """ Calculate normalized frequency for a given temperature and frequency.

    :param temp: Temperature.
    :param omega: Frequency.
    :returns Distribution which is needed for determining the heat capacity.
    """
    return inv_cm_2_ev / k_2_ev * omega / temp if temp > eps_min else eps_max

def vibrational_expression(omega: float) -> float:
    """ Calculate distribution which is later multiplied by the phonon DOS to determine the vibrational free energy.

    :param omega: Frequency.
    :returns Distribution which is needed for determining the vibrational free energy.
    """
    omega_dist = 1 / 2
    if omega > eps_min:
        omega_dist = omega / 2
        if omega < eps_max:
            omega_dist = omega / 2 + np.log(1 - np.exp(-omega))

    return omega_dist

def thermodynamic_properties(
        temp: float, omega_data: np.ndarray, dos_data: np.ndarray, energy_unit: str
) -> Tuple[float, float, float, float, float, float]:
    """ Calculate thermodynamic properties at a given temperature.

    :param temp: Temperature.
    :param omega_data: Frequency values.
    :param dos_data: DOS values.
    :param energy_unit: Energy unit.
    :returns: Vibrational free energy, vibrational internal energy, entropic contribution to the vibrational free
    energy, vibrational entropy, heat capacity, zero-point energy
    """
    conversion_factor = ev_2_ha if energy_unit == 'Ha' else 1
    free_energy = heat_capacity = zero_point_energy = internal_energy = entropy = 0

    for i in range(len(dos_data) - 1):
        delta_omega = omega_data[i + 1] - omega_data[i]
        x = normalized_frequency(temp, omega_data[i])
        free_energy += delta_omega * dos_data[i] * vibrational_expression(x) * temp * k_2_ev * conversion_factor
        heat_capacity += delta_omega * dos_data[i] * heat_capacity_distribution(x)
        zero_point_energy += delta_omega * dos_data[i] * x * temp * k_2_ev * conversion_factor / 2
        internal_energy += delta_omega * dos_data[i] * x / np.tanh(x / 2) * temp * k_2_ev * conversion_factor / 2
        entropy += delta_omega * dos_data[i] * (x / np.tanh(x / 2) / 2 - vibrational_expression(x)) \
                   * conversion_factor * k_2_ev * 1000

    if temp < 1.e-6: entropy = 0
    ts_vib = internal_energy - free_energy

    return free_energy, internal_energy, ts_vib, entropy, heat_capacity, zero_point_energy


def main() -> None:
    parser = ArgumentParser(description="""Calculate thermodynamic properties from the phonon DOS.
                                        Temperatures should be given in Kelvin""")

    parser.add_argument("t_min",
                        type=float,
                        nargs=1,
                        help="minimum temperature")

    parser.add_argument("t_max",
                        type=float,
                        nargs=1,
                        help="maximum temperature")

    parser.add_argument("n_steps",
                        type=int,
                        nargs=1,
                        help="number of steps in temperature range")

    parser.add_argument("energy_unit",
                        type=str,
                        nargs=1,
                        default=["eV"],
                        help="energy unit")

    args = parser.parse_args()


    t_step = (args.t_max[0] - args.t_min[0]) / args.n_steps[0]
    temperature_range = [max(args.t_min[0] + i * t_step, 1e-10) for i in range(args.n_steps[0] + 1)]

    phonon_dos_data = np.genfromtxt("PHDOS.OUT")
    omega_data = phonon_dos_data[:, 0] * ha_2_inv_cm
    dos_data = phonon_dos_data[:, 1] / ha_2_inv_cm

    file_names = ['F_vib', 'U_vib', 'TS_vib', 'S_vib', 'C_v']
    files = [open(file_name, "w") for file_name in file_names]

    for temp in temperature_range:
        properties = thermodynamic_properties(temp, omega_data, dos_data, args.energy_unit[0])
        for index, value in enumerate(properties[:5]):
            files[index].write(f'{temp:.8e} {value:.8e}\n')

        if temp == temperature_range[0]:
            print(f"\nZero-point energy is {properties[5]:.4e} [{args.energy_unit[0]}]\n")

    for file in files:
        file.close()


if __name__ == "__main__":
    main()