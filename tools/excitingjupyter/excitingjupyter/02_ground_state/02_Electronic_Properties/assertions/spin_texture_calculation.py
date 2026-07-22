from os.path import dirname

import numpy as np

from excitingtools.exciting_dict_parsers.groundstate_parser import parse_info_out

TUTORIAL_RUNDIR = "../run_tutorial_spin_texture"


def test_groundstate(converged_results: dict):
    """Test results for groundstate calculation in main output file INFO.OUT.
    """

    total_energy = -4205.16171282
    assert np.isclose(converged_results['Total energy'], total_energy), \
        f"Incorrect value for total energy. {total_energy} != {converged_results['Total energy']}"

    fermi_energy = 0.18407192
    assert np.isclose(converged_results['Fermi energy'], fermi_energy), \
        f"Incorrect value for fermi energy. {fermi_energy} != {converged_results['Fermi energy']}"

    kinetic_energy = 4294.91833245
    assert np.isclose(converged_results['Kinetic energy'], kinetic_energy), \
        f"Incorrect value for kinetic energy. {kinetic_energy} != {converged_results['Kinetic energy']}"

    coulomb_energy = -8340.45886496
    assert np.isclose(converged_results['Coulomb energy'], coulomb_energy), \
        f"Incorrect value for coulomb energy. {coulomb_energy} != {converged_results['Coulomb energy']}"

    exchange_energy = -156.49570126
    assert np.isclose(converged_results['Exchange energy'], exchange_energy), \
        f"Incorrect value for exchange energy. {exchange_energy} != {converged_results['Exchange energy']}"

    correlation_energy = -3.12547905
    assert np.isclose(converged_results['Correlation energy'], correlation_energy), \
        f"Incorrect value for correlation energy. {correlation_energy} != {converged_results['Correlation energy']}"

    fermi_energy_dos = 0.00000003
    assert np.isclose(converged_results['DOS at Fermi energy (states/Ha/cell)'], fermi_energy_dos), \
        f"Incorrect value for DOS at Fermi energy (states/Ha/cell). {fermi_energy_dos} != {converged_results['DOS at Fermi energy (states/Ha/cell)']}"

    charge_core_electrons = 36.00000000
    assert np.isclose(converged_results['core'], charge_core_electrons), \
        f"Incorrect value for charge of core electrons. {charge_core_electrons} != {converged_results['core']}"

    core_leakage_charge = 0.00765189
    assert np.isclose(converged_results['core leakage'], core_leakage_charge), \
        f"Incorrect value for core leakage charge. {core_leakage_charge} != {converged_results['core leakage']}"

    charge_valance_electrons = 28.00000000
    assert np.isclose(converged_results['valence'], charge_valance_electrons), \
        f"Incorrect value for charge of valence electrons. {charge_valance_electrons} != {converged_results['valence']}"

    interstitial_region_charge = np.array([-0.0, 0.0, -0.0])
    assert np.allclose(converged_results['interstitial'], interstitial_region_charge), \
        f"Incorrect value for charge in interstitial region. {interstitial_region_charge} != {converged_results['interstitial']}"

    muffin_tins_charge = 59.15523680
    assert np.isclose(converged_results['total charge in muffin-tins'], muffin_tins_charge), \
        f"Incorrect value for total charge in muffin-tins. {muffin_tins_charge} != {converged_results['total charge in muffin-tins']}"

    total_charge = 64.00000000
    assert np.isclose(converged_results['total charge'], total_charge), \
        f"Incorrect value for total charge. {total_charge} != {converged_results['total charge']}"

    estimated_gap = 0.01393823
    assert np.isclose(converged_results['Estimated fundamental gap'], estimated_gap, atol=1e-5), \
        f"Incorrect value for Estimated fundamental gap. {estimated_gap} != {converged_results['Estimated fundamental gap']}"


def main():
    results = parse_info_out(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/INFO.OUT")
    max_scf = max([int(i) for i in results['scl'].keys()])
    assert max_scf <= 20, "Expect max 20 SCF iterations to converge"
    converged_results = results['scl'][str(max_scf)]
    test_groundstate(converged_results)

if __name__=="__main__":
    main()
