from os.path import dirname

import numpy as np

TUTORIAL_RUNDIR = "../run_tutorial_van_der_waals_corrections"


def test_van_der_waals():
    energy_vs_strain_DFTD2_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/DFTD2/energy-vs-strain")
    energy_vs_strain_DFTD2_reference = np.array([[5.00000000, -152.41224082],
                                                [5.33333333, -152.42211038],
                                                [5.66666667, -152.42633611],
                                                [6.00000000, -152.42761455],
                                                [6.33333333, -152.42750741],
                                                [6.66666667, -152.42684059],
                                                [7.00000000, -152.42597225],
                                                [7.33333333, -152.42512810],
                                                [7.66666667, -152.42436899],
                                                [8.00000000, -152.42367094],
                                                [20.00000000, -152.41932763]])

    assert np.allclose(energy_vs_strain_DFTD2_results, energy_vs_strain_DFTD2_reference,  atol=1.0e-8),(
                       "DFTD2 energy vs strain values not equivalent to reference calculation")

    energy_vs_strain_PBE_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/PBE-uncorrected/energy-vs-strain")
    energy_vs_strain_PBE_reference = np.array([[5.00000000, -152.38457377],
                                              [5.33333333, -152.39636527],
                                              [5.66666667, -152.40320916],
                                              [6.00000000, -152.40713649],
                                              [6.33333333, -152.40933239],
                                              [6.66666667, -152.41052651],
                                              [7.00000000, -152.41113087],
                                              [7.33333333, -152.41145310],
                                              [7.66666667, -152.41162457],
                                              [8.00000000, -152.41167557],
                                              [20.00000000, -152.41126582]])

    assert np.allclose(energy_vs_strain_PBE_results, energy_vs_strain_PBE_reference, atol=1.0e-8),(
                       "PBE energy vs strain values not equivalent to reference calculation")

    energy_vs_strain_TSvdW_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/TSvdW/energy-vs-strain")
    energy_vs_strain_TSvdW_reference = np.array([[5.00000000, -152.40996062],
                                                [5.33333333, -152.42048122],
                                                [5.66666667, -152.42649168],
                                                [6.00000000, -152.42954783],
                                                [6.33333333, -152.43030636],
                                                [6.66666667, -152.42960948],
                                                [7.00000000, -152.42827166],
                                                [7.33333333, -152.42686796],
                                                [7.66666667, -152.42560714],
                                                [8.00000000, -152.42449556],
                                                [20.00000000, -152.41818978]])

    assert np.allclose(energy_vs_strain_TSvdW_results, energy_vs_strain_TSvdW_reference, atol=1.0e-8),(
                       "TSvdW energy vs strain values not equivalent to reference calculation")


def main():
    test_van_der_waals()

if __name__=="__main__":
    main()
