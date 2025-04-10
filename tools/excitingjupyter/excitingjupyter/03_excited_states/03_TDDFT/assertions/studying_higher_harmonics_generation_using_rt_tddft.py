from os.path import dirname

import numpy as np

TUTORIAL_RUNDIR = "../run_MoS2_hhg/ampl_1"
REFERENCE_DIR = "reference_studying_higher_harmonics_generation_using_rt_tddft"


def test_higher_harmonic_generation_rt_tddft():
    pvec_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/PVEC.OUT")
    pvec_reference = np.loadtxt(f"{dirname(__file__)}/{REFERENCE_DIR}/PVEC.OUT")

    assert np.allclose(pvec_results, pvec_reference,  atol=1.0e-7),(
                       "Polarization vector results not equivalent to reference calculation")


def main():
    test_higher_harmonic_generation_rt_tddft()

if __name__=="__main__":
    main()
