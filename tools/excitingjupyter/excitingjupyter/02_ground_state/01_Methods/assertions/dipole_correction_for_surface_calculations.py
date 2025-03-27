from os.path import dirname

import numpy as np

from excitingtools.exciting_dict_parsers.groundstate_parser import parse_info_out

TUTORIAL_RUNDIR = "../run_tutorial_dipole_correction/DFT-dc"


def test_groundstate(converged_results: dict):
    """Test results for groundstate calculation in main output file INFO.OUT.
    """

    total_energy = -175.49540334
    assert np.isclose(converged_results['Total energy'], total_energy), \
        f"Incorrect value for total energy. Expect {total_energy}"

    fermi_energy = -0.07775235
    assert np.isclose(converged_results['Fermi energy'], fermi_energy), \
        f"Incorrect value for fermi energy. Expect {fermi_energy}"

    kinetic_energy = 174.46531741
    assert np.isclose(converged_results['Kinetic energy'], kinetic_energy), \
        f"Incorrect value for kinetic energy. Expect {kinetic_energy}"

    coulomb_energy = -329.77860949
    assert np.isclose(converged_results['Coulomb energy'], coulomb_energy), \
        f"Incorrect value for coulomb energy. Expect {coulomb_energy}"

    exchange_energy = -18.81370447
    assert np.isclose(converged_results['Exchange energy'], exchange_energy), \
        f"Incorrect value for exchange energy. Expect {exchange_energy}"

    correlation_energy = -1.49667048
    assert np.isclose(converged_results['Correlation energy'], correlation_energy), \
        f"Incorrect value for correlation energy. Expect {correlation_energy}"

    fermi_energy_dos = 0.00000000
    assert np.isclose(converged_results['DOS at Fermi energy (states/Ha/cell)'], fermi_energy_dos), \
        f"Incorrect value for DOS at Fermi energy (states/Ha/cell). Expect {fermi_energy_dos}"

    charge_core_electrons = 6.00000000
    assert np.isclose(converged_results['core'], charge_core_electrons), \
        f"Incorrect value for charge of core electrons. Expect {charge_core_electrons}"

    core_leakage_charge = 0.00090344
    assert np.isclose(converged_results['core leakage'], core_leakage_charge), \
        f"Incorrect value for core leakage charge. Expect {core_leakage_charge}"

    charge_valance_electrons = 16.00000000
    assert np.isclose(converged_results['valence'], charge_valance_electrons), \
        f"Incorrect value for charge of valence electrons. Expect {charge_valance_electrons}"

    interstitial_region_charge = 7.40332722
    assert np.isclose(converged_results['interstitial'], interstitial_region_charge), \
        f"Incorrect value for charge in interstitial region. Expect {interstitial_region_charge}"

    muffin_tins_charge = 14.59667278
    assert np.isclose(converged_results['total charge in muffin-tins'], muffin_tins_charge), \
        f"Incorrect value for total charge in muffin-tins. Expect {muffin_tins_charge}"

    total_charge = 22.00000000
    assert np.isclose(converged_results['total charge'], total_charge), \
        f"Incorrect value for total charge. Expect {total_charge}"

    estimated_gap = 0.10519764
    assert np.isclose(converged_results['Estimated fundamental gap'], estimated_gap), \
        f"Incorrect value for Estimated fundamental gap. Expect {estimated_gap}"


def main():
    results = parse_info_out(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/INFO.OUT")
    max_scf = max([int(i) for i in results['scl'].keys()])
    assert max_scf <= 20, "Expect max 20 SCF iterations to converge"
    converged_results = results['scl'][str(max_scf)]
    test_groundstate(converged_results)

if __name__=="__main__":
    main()
