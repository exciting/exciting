from os.path import dirname

import numpy as np

TUTORIAL_RUNDIR = "../run_BN_Ehrenfest/ehrenfest"
REFERENCE_DIR = "reference_real_time_tddft_combined_with_molecular_dynamics"


def test_rt_tddft_md():
    atom_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/ATOM_0001.OUT")
    atom_reference = np.loadtxt(f"{dirname(__file__)}/{REFERENCE_DIR}/ATOM_0001.OUT")

    assert np.allclose(atom_results, atom_reference,  atol=1.0e-7),(
                       "Ehrenfest Dynamics results not equivalent to reference calculation")


def main():
    test_rt_tddft_md()

if __name__=="__main__":
    main()
