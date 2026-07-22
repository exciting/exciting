from os.path import dirname

import numpy as np

from excitingtools.exciting_dict_parsers.groundstate_parser import parse_info_out

TUTORIAL_RUNDIR = "../run_graphene"
REFERENCE_DIR = "reference_graphene"


def test_groundstate(converged_results: dict):
    """Test results for groundstate calculation in main output file INFO.OUT.
    """

    total_energy = -76.11387789
    assert np.isclose(converged_results['Total energy'], total_energy), \
        f"Incorrect value for total energy. Expect {total_energy}"

    fermi_energy = -0.06783556
    assert np.isclose(converged_results['Fermi energy'], fermi_energy), \
        f"Incorrect value for fermi energy. Expect {fermi_energy}"

    kinetic_energy = 75.81190852
    assert np.isclose(converged_results['Kinetic energy'], kinetic_energy), \
        f"Incorrect value for kinetic energy. Expect {kinetic_energy}"

    coulomb_energy = -141.13032000
    assert np.isclose(converged_results['Coulomb energy'], coulomb_energy), \
        f"Incorrect value for coulomb energy. Expect {coulomb_energy}"

    exchange_energy = -10.37420659
    assert np.isclose(converged_results['Exchange energy'], exchange_energy), \
        f"Incorrect value for exchange energy. Expect {exchange_energy}"

    correlation_energy = -0.42125982
    assert np.isclose(converged_results['Correlation energy'], correlation_energy), \
        f"Incorrect value for correlation energy. Expect {correlation_energy}"

    fermi_energy_dos = 0.00002749
    assert np.isclose(converged_results['DOS at Fermi energy (states/Ha/cell)'], fermi_energy_dos), \
        f"Incorrect value for DOS at Fermi energy (states/Ha/cell). Expect {fermi_energy_dos}"

    charge_core_electrons = 4.00000000
    assert np.isclose(converged_results['core'], charge_core_electrons), \
        f"Incorrect value for charge of core electrons. Expect {charge_core_electrons}"

    core_leakage_charge = 0.00129555
    assert np.isclose(converged_results['core leakage'], core_leakage_charge), \
        f"Incorrect value for core leakage charge. Expect {core_leakage_charge}"

    charge_valance_electrons = 8.00000000
    assert np.isclose(converged_results['valence'], charge_valance_electrons), \
        f"Incorrect value for charge of valence electrons. Expect {charge_valance_electrons}"

    interstitial_region_charge = 4.85400965
    assert np.isclose(converged_results['interstitial'], interstitial_region_charge), \
        f"Incorrect value for charge in interstitial region. Expect {interstitial_region_charge}"

    muffin_tins_charge = 7.14599035
    assert np.isclose(converged_results['total charge in muffin-tins'], muffin_tins_charge), \
        f"Incorrect value for total charge in muffin-tins. Expect {muffin_tins_charge}"

    total_charge = 12.00000000
    assert np.isclose(converged_results['total charge'], total_charge), \
        f"Incorrect value for total charge. Expect {total_charge}"

    estimated_gap = 0.07363036
    assert np.isclose(converged_results['Estimated fundamental gap'], estimated_gap), \
        f"Incorrect value for Estimated fundamental gap. Expect {estimated_gap}"

def test_energy_vs_strain():
    energy_vs_strain_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/opt/opt-run/energy-vs-strain")
    energy_vs_strain_reference = np.array([[-0.10000000, -76.03156387],
                                        [-0.09000000, -76.05036272],
                                        [-0.08000000, -76.06602406],
                                        [-0.07000000, -76.07894462],
                                        [-0.06000000, -76.08890082],
                                        [-0.05000000, -76.09744293],
                                        [-0.04000000, -76.10354890],
                                        [-0.03000000, -76.10908624],
                                        [-0.02000000, -76.11277032],
                                        [-0.01000000, -76.11356163],
                                        [ 0.00000000, -76.11387789],
                                        [ 0.01000000, -76.11410361],
                                        [ 0.02000000, -76.11163222],
                                        [ 0.03000000, -76.10909896],
                                        [ 0.04000000, -76.10729404],
                                        [ 0.05000000, -76.10287673],
                                        [ 0.06000000, -76.09885160],
                                        [ 0.07000000, -76.09398243],
                                        [ 0.08000000, -76.08854410],
                                        [ 0.09000000, -76.08276160],
                                        [ 0.10000000, -76.07716831]])

    assert np.allclose(energy_vs_strain_results, energy_vs_strain_reference, atol=1.0e-7), (
        "energy vs strain values not equivalent to reference calculation")

def test_energy_vs_displacement():
    energy_vs_displacement_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/planar/str-0.10/energy-vs-displacement")
    energy_vs_displacement_reference = np.array([[-0.20000000, -75.7434488215],
                                                [-0.18000000, -75.7500826157],
                                                [-0.16000000, -75.7638450774],
                                                [-0.14000000, -75.7824895650],
                                                [-0.12000000, -75.8024964249],
                                                [-0.10000000, -75.8228567348],
                                                [-0.08000000, -75.8623631368],
                                                [-0.06000000, -75.9239703391],
                                                [-0.04000000, -76.0221798077],
                                                [-0.02000000, -76.0961349318],
                                                [ 0.00000100, -76.1208558516],
                                                [ 0.02000000, -76.0961349318],
                                                [ 0.04000000, -76.0221798077],
                                                [ 0.06000000, -75.9239703391],
                                                [ 0.08000000, -75.8623631368],
                                                [ 0.10000000, -75.8228567348],
                                                [ 0.12000000, -75.8024964249],
                                                [ 0.14000000, -75.7824895650],
                                                [ 0.16000000, -75.7638450774],
                                                [ 0.18000000, -75.7500826157],
                                                [ 0.20000000, -75.7434488215]])

    assert np.allclose(energy_vs_displacement_results, energy_vs_displacement_reference, atol=1.0e-7), (
        "energy vs displacement values not equivalent to reference calculation")

def test_loss():
    loss_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/loss-function/loss-k09")
    loss_reference = np.loadtxt(f"{dirname(__file__)}/{REFERENCE_DIR}/loss-k09")

    assert np.allclose(loss_results, loss_reference,  atol=1.0e-7),(
                       "Loss function results not equivalent to reference calculation")

def main():
    results = parse_info_out(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/test-1/INFO.OUT")
    max_scf = max([int(i) for i in results['scl'].keys()])
    assert max_scf <= 25, "Expect max 25 SCF iterations to converge"
    converged_results = results['scl'][str(max_scf)]
    test_groundstate(converged_results)

    test_energy_vs_strain()
    test_energy_vs_displacement()
    test_loss()

if __name__=="__main__":
    main()
