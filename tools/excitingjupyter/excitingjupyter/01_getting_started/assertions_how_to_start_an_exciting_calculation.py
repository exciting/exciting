import numpy as np
from os.path import dirname

from excitingtools import parser_chooser

TUTORIAL_HOW_TO_START_AN_EXCITING_CALCULATION_RUNDIR = "run_tutorial_start_exciting_calculation"

def test_tutorial1(converged_results):
    """Automatically test results of 01_getting_started notebook, diamond bulk calculation.
    """
    assert np.isclose(converged_results['Total energy'], -75.88903685), (
        f"Total energy in Ha is {converged_results['Total energy']}")

    assert np.isclose(converged_results['Fermi energy'], 0.50878171), (
        f"Fermi energy in Ha is {converged_results['Fermi energy']}")

    assert np.isclose(converged_results['Kinetic energy'], 75.54873237), (
        f"Kinetic energy in Ha is {converged_results['Kinetic energy']}")

    assert np.isclose(converged_results['Coulomb energy'], -140.90113491), (
        f"Coulomb energy in Ha is {converged_results['Coulomb energy']}")

    assert np.isclose(converged_results['Exchange energy'], -10.03412627), (
        f"Exchange energy in Ha is {converged_results['Exchange energy']}")

    assert np.isclose(converged_results['Correlation energy'], -0.50250803), (
        f"Correlation energy in Ha is {converged_results['Correlation energy']}")

    assert np.isclose(converged_results['DOS at Fermi energy (states/Ha/cell)'], 0.0), (
        f"DOS at Fermi energy in states/Ha/cell is "
        f"{converged_results['DOS at Fermi energy (states/Ha/cell)']}")

    assert np.isclose(converged_results['core'], 4.0), (
        f"Charge of core electrons is {converged_results['core']}")

    assert np.isclose(converged_results['core leakage'], 0.00016223), (
        f"Core leakage charge is {converged_results['core leakage']}")

    assert np.isclose(converged_results['valence'], 8.0), (
        f"Charge of valence electrons is {converged_results['valence']}")

    assert np.isclose(converged_results['interstitial'], 3.09358006), (
        f"Charge in the interstitial region is {converged_results['interstitial']}")

    assert np.isclose(converged_results['atom     1     C'], 4.45320997), (
        f"Charge in muffin-tin spheres: 1st stom (C) is {converged_results['atom     1     C']}")

    assert np.isclose(converged_results['atom     2     C'], 4.45320997), (
        f"Charge in muffin-tin spheres: 2nd stom (C) is {converged_results['atom     2     C']}")

    assert np.isclose(converged_results['total charge in muffin-tins'], 8.90641994), (
        f"Total charge in muffin-tins is {converged_results['total charge in muffin-tins']}")

    assert np.isclose(converged_results['total charge'], 12.0), (
        f"Total charge is {converged_results['total charge']}")

    assert np.isclose(converged_results['Estimated fundamental gap'], 0.15982903), (
        f"Estimated fundamental gap in Ha is {converged_results['Estimated fundamental gap']}")

def main():
    results = parser_chooser(f"{dirname(__file__)}/{TUTORIAL_HOW_TO_START_AN_EXCITING_CALCULATION_RUNDIR}/INFO.OUT")
    max_scf = max([int(i) for i in results['scl'].keys()])
    assert max_scf <= 13, "Expect max 13 SCF iterations to converge"
    converged_results = results['scl'][str(max_scf)]

    test_tutorial1(converged_results)

if __name__=="__main__":
    main()