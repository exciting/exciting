from os.path import dirname

import numpy as np

TUTORIAL_RUNDIR = "../run_BE_Si"
REFERENCE_DIR = "reference_transport_properties_using_the_boltzmann_equation"


def test_boltzmann():
    seebeck_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/SEEBECK_11.OUT")
    seebeck_reference = np.loadtxt(f"{dirname(__file__)}/{REFERENCE_DIR}/SEEBECK_11.OUT")

    assert np.allclose(seebeck_results, seebeck_reference,  atol=1.0e-7),(
                       "Seebeck coefficient values not equivalent to reference calculation")

def main():
    test_boltzmann()

if __name__=="__main__":
    main()