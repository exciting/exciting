from os.path import dirname

import numpy as np

TUTORIAL_RUNDIR = "../moke/k444"
REFERENCE_DIR = "reference_magneto_optical_kerr_effect"


def test_moke():
    moke_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/MOKE_NLF_FXCRPA_QMT001.OUT")
    moke_reference = np.loadtxt(f"{dirname(__file__)}/{REFERENCE_DIR}/MOKE_NLF_FXCRPA_QMT001.OUT")

    assert np.allclose(moke_results, moke_reference,  atol=1.0e-7),(
                       "MOKE values not equivalent to reference calculation")

def main():
    test_moke()

if __name__=="__main__":
    main()