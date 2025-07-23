from os.path import dirname

import numpy as np

TUTORIAL_RUNDIR = "../run_MoS2_pump_probe"
REFERENCE_DIR = "reference_simulating_pump_probe_spectroscopy_with_rt_tddft"


def test_pump_probe():
    jind_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/DIFF.OUT")
    jind_reference = np.loadtxt(f"{dirname(__file__)}/{REFERENCE_DIR}/JIND_DIFF.OUT")

    assert np.allclose(jind_results, jind_reference,  atol=1.0e-7),(
                       "Current density as the difference between the values after the pump and the probe and the \
                        pump pulse  not equivalent to reference calculation")


def main():
    test_pump_probe()

if __name__=="__main__":
    main()
