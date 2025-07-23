from os.path import dirname

import numpy as np

TUTORIAL_RUNDIR = "../shg"
REFERENCE_DIR = "reference_second_harmonic_generation"


def test_shg():
    chi_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/CHI_123.OUT")
    chi_reference = np.loadtxt(f"{dirname(__file__)}/{REFERENCE_DIR}/CHI_123.OUT")

    assert np.allclose(chi_results, chi_reference,  atol=1.0e-6),(
                       "Susceptibility tensor values not equivalent to reference calculation")

def main():
    test_shg()

if __name__=="__main__":
    main()