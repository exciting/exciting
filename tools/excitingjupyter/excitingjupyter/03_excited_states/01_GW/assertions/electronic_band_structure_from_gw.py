from os.path import dirname

import numpy as np
from excitingtools import parse

TUTORIAL_GW_RUNDIR = "run_Si_GW"
REFERENCE_DIR = "reference_electronic_band_structure_from_gw"


def test_groundstate(converged_results: dict):
    """Test results for ground-state calculation in main output file INFO.OUT.
    """

    total_energy = -578.05734534
    assert np.isclose(converged_results['Total energy'], total_energy), \
        f"Incorrect value for total energy. Expect {total_energy}"

    fermi_energy = 0.20328046
    assert np.isclose(converged_results['Fermi energy'], fermi_energy), \
        f"Incorrect value for fermi energy. Expect {fermi_energy}"

    kinetic_energy = 579.09821116
    assert np.isclose(converged_results['Kinetic energy'], kinetic_energy), \
        f"Incorrect value for kinetic energy. Expect {kinetic_energy}"

    coulomb_energy = -1117.45646299
    assert np.isclose(converged_results['Coulomb energy'], coulomb_energy), \
        f"Incorrect value for coulomb energy. Expect {coulomb_energy}"

    exchange_energy = -37.55285876
    assert np.isclose(converged_results['Exchange energy'], exchange_energy), \
        f"Incorrect value for exchange energy. Expect {exchange_energy}"

    correlation_energy = -2.14623474
    assert np.isclose(converged_results['Correlation energy'], correlation_energy), \
        f"Incorrect value for correlation energy. Expect {correlation_energy}"

    fermi_energy_dos = 0.00000000
    assert np.isclose(converged_results['DOS at Fermi energy (states/Ha/cell)'], fermi_energy_dos), \
        f"Incorrect value for DOS at Fermi energy (states/Ha/cell). Expect {fermi_energy_dos}"

    charge_core_electrons = 20.00000000
    assert np.isclose(converged_results['core'], charge_core_electrons), \
        f"Incorrect value for charge of core electrons. Expect {charge_core_electrons}"

    core_leakage_charge = 0.00461386
    assert np.isclose(converged_results['core leakage'], core_leakage_charge), \
        f"Incorrect value for core leakage charge. Expect {core_leakage_charge}"

    charge_valance_electrons = 8.00000000
    assert np.isclose(converged_results['valence'], charge_valance_electrons), \
        f"Incorrect value for charge of valence electrons. Expect {charge_valance_electrons}"

    interstitial_region_charge = 3.87731748
    assert np.isclose(converged_results['interstitial'], interstitial_region_charge), \
        f"Incorrect value for charge in interstitial region. Expect {interstitial_region_charge}"

    muffin_tins_charge = 24.12268252
    assert np.isclose(converged_results['total charge in muffin-tins'], muffin_tins_charge), \
        f"Incorrect value for total charge in muffin-tins. Expect {muffin_tins_charge}"

    total_charge = 28.00000000
    assert np.isclose(converged_results['total charge'], total_charge), \
        f"Incorrect value for total charge. Expect {total_charge}"

    estimated_gap = 0.02094574
    assert np.isclose(converged_results['Estimated fundamental gap'], estimated_gap), \
        f"Incorrect value for Estimated fundamental gap. Expect {estimated_gap}"


def test_band_ks(file_name: str):
    """Test results of Kohn-Sham (KS) electronic band structure calculations.
    """

    reference_band_ks = np.genfromtxt(f"{dirname(__file__)}/{REFERENCE_DIR}/{file_name}")
    band_results = np.genfromtxt(f"{dirname(__file__)}/../{TUTORIAL_GW_RUNDIR}/{file_name}")

    assert np.allclose(band_results[:, 0], reference_band_ks[:, 0]), \
        f"K-grid of the Kohn-Sham (KS) electronic band-structure calculation along high-symmetry lines in the" \
        f" Brillouin zone not equivalent to reference calculation."

    assert np.allclose(band_results[:, 1], reference_band_ks[:, 1]), \
        f"Energy values of the Kohn-Sham (KS) electronic band-structure calculation not equivalent to reference" \
        f" calculation."


def test_dos_ks(file_name: str):
    """Test results of Kohn-Sham (KS) density of states calculations.
    """

    reference_dos_ks = np.genfromtxt(f"{dirname(__file__)}/{REFERENCE_DIR}/{file_name}")
    dos_results = np.genfromtxt(f"{dirname(__file__)}/../{TUTORIAL_GW_RUNDIR}/{file_name}")

    assert np.allclose(dos_results[:, 0], reference_dos_ks[:, 0]), \
        f"Energy values the Kohn-Sham (KS) density of states calculation not equivalent to reference calculation."

    assert np.allclose(dos_results[:, 1], reference_dos_ks[:, 1]), \
        f"Kohn-Sham (KS) density of states(states/eV/unit cell) values not equivalent to reference calculation."


# GW will give results that vary on the order ~ 1 meV per QP eigenvalue if MKL thread settings are not consistent,
# and the calculation does not start from STATE.OUT. While this is physically small, it prevents testing the results

def main():
    results = parse(f"{dirname(__file__)}/../{TUTORIAL_GW_RUNDIR}/INFO.OUT")
    max_scf = max([int(i) for i in results['scl'].keys()])
    assert max_scf <= 12, "Expect max 12 SCF iterations to converge"
    converged_results = results['scl'][str(max_scf)]
    test_groundstate(converged_results)

    test_band_ks("BAND.OUT")

    test_dos_ks("TDOS.OUT")

if __name__=="__main__":
    main()
