from os.path import dirname

import numpy as np

TUTORIAL_RUNDIR = "../run_cBN_Borncharges/dfpt"


def test_born_effective_charges():
    zstar_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/ZSTAR.OUT")
    zstar_reference = np.array([[ 2.0052078e+00,  0.0000000e+00,  0.0000000e+00],
                                [ 0.0000000e+00,  2.0052078e+00,  0.0000000e+00],
                                [ 0.0000000e+00,  0.0000000e+00,  2.0052078e+00],
                                [-2.0052078e+00,  0.0000000e+00,  0.0000000e+00],
                                [ 0.0000000e+00, -2.0052078e+00,  0.0000000e+00],
                                [ 0.0000000e+00,  0.0000000e+00, -2.0052078e+00],
                                [ 2.6674000e-04,  0.0000000e+00,  0.0000000e+00],
                                [ 0.0000000e+00,  2.6674000e-04,  0.0000000e+00],
                                [ 0.0000000e+00,  0.0000000e+00,  2.6674000e-04]])

    assert np.allclose(zstar_results, zstar_reference,  atol=1.0e-8),(
                       "energy vs volume values not equivalent to reference calculation")



def main():
    test_born_effective_charges()

if __name__=="__main__":
    main()
