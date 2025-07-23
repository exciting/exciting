from os.path import dirname

import numpy as np

from excitingtools.exciting_dict_parsers.groundstate_parser import parse_info_out

TUTORIAL_RUNDIR = "../CO2/HOMO_2"


def test_groundstate(converged_results: dict):
    """Test results for groundstate calculation in main output file INFO.OUT.
    """

    total_energy = -186.67876871
    assert np.isclose(converged_results['Total energy'], total_energy), \
        f"Incorrect value for total energy. Expect {total_energy}"

    fermi_energy = 0.08314148
    assert np.isclose(converged_results['Fermi energy'], fermi_energy), \
        f"Incorrect value for fermi energy. Expect {fermi_energy}"

    fermi_energy_dos = 0.00000000
    assert np.isclose(converged_results['DOS at Fermi energy (states/Ha/cell)'], fermi_energy_dos), \
        f"Incorrect value for DOS at Fermi energy (states/Ha/cell). Expect {fermi_energy_dos}"

    estimated_gap = 0.14821665
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
