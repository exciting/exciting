from os.path import dirname

import numpy as np

from excitingtools.exciting_dict_parsers.groundstate_parser import parse_info_out

TUTORIAL_RUNDIR = "../run_diamond_rt"

def test_real_time_tddft(converged_results):
    """Test results for RT-TDDFT calculation.
    """
    total_energy = -75.61687858
    assert np.isclose(converged_results['Total energy'], total_energy), \
        f"Incorrect value for total energy. Expect {total_energy}"

    fermi_energy = 0.49637627
    assert np.isclose(converged_results['Fermi energy'], fermi_energy), \
        f"Incorrect value for fermi energy. Expect {fermi_energy}"

    kinetic_energy = 75.03596256
    assert np.isclose(converged_results['Kinetic energy'], kinetic_energy), \
        f"Incorrect value for kinetic energy. Expect {kinetic_energy}"

    coulomb_energy = -140.41345770
    assert np.isclose(converged_results['Coulomb energy'], coulomb_energy), \
        f"Incorrect value for coulomb energy. Expect {coulomb_energy}"

    exchange_energy = -9.43820947
    assert np.isclose(converged_results['Exchange energy'], exchange_energy), \
        f"Incorrect value for exchange energy. Expect {exchange_energy}"

    correlation_energy = -0.80117397
    assert np.isclose(converged_results['Correlation energy'], correlation_energy), \
        f"Incorrect value for correlation energy. Expect {correlation_energy}"

    fermi_energy_dos = 0.00000000
    assert np.isclose(converged_results['DOS at Fermi energy (states/Ha/cell)'], fermi_energy_dos), \
        f"Incorrect value for DOS at Fermi energy (states/Ha/cell). Expect {fermi_energy_dos}"

    charge_core_electrons = 4.00000000
    assert np.isclose(converged_results['core'], charge_core_electrons), \
        f"Incorrect value for charge of core electrons. Expect {charge_core_electrons}"

    core_leakage_charge = 0.00025692
    assert np.isclose(converged_results['core leakage'], core_leakage_charge), \
        f"Incorrect value for core leakage charge. Expect {core_leakage_charge}"

    charge_valance_electrons = 8.00000000
    assert np.isclose(converged_results['valence'], charge_valance_electrons), \
        f"Incorrect value for charge of valence electrons. Expect {charge_valance_electrons}"

    interstitial_region_charge = 3.49527867
    assert np.isclose(converged_results['interstitial'], interstitial_region_charge), \
        f"Incorrect value for charge in interstitial region. Expect {interstitial_region_charge}"

    muffin_tins_charge = 8.50472133
    assert np.isclose(converged_results['total charge in muffin-tins'], muffin_tins_charge), \
        f"Incorrect value for total charge in muffin-tins. Expect {muffin_tins_charge}"

    total_charge = 12.00000000
    assert np.isclose(converged_results['total charge'], total_charge), \
        f"Incorrect value for total charge. Expect {total_charge}"

    estimated_gap = 0.16259259
    assert np.isclose(converged_results['Estimated fundamental gap'], estimated_gap), \
        f"Incorrect value for Estimated fundamental gap. Expect {estimated_gap}"

def main():
    results = parse_info_out(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/INFO_RTTDDFT_GND.OUT")
    max_scf = max([int(i) for i in results['scl'].keys()])
    converged_results = results['scl'][str(max_scf)]

    test_real_time_tddft(converged_results)

if __name__=="__main__":
    main()