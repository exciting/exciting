from os.path import dirname

import numpy as np

TUTORIAL_RUNDIR = "../silver-volume-optimization/vol-opt-1"


def test_volume_optimization_cubic_systems():
    energy_vs_volume_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/energy-vs-volume")
    energy_vs_volume_reference = np.array([[98.9647988854, -5314.77888971],
                                           [102.1230153745, -5314.78069604],
                                           [105.3477182017, -5314.78200268],
                                           [108.639599933, -5314.78329031],
                                           [111.9993531345, -5314.78341370],
                                           [115.4280166556, -5314.78296499],
                                           [118.9252442122, -5314.78284893],
                                           [122.4927672204, -5314.78187065],
                                           [126.1309319629, -5314.78079209],
                                           [129.8404310056, -5314.77916096],
                                           [133.6219569147, -5314.77829916]])

    assert np.allclose(energy_vs_volume_results, energy_vs_volume_reference,  atol=1.0e-8),(
                       "energy vs volume values not equivalent to reference calculation")



def main():
    test_volume_optimization_cubic_systems()

if __name__=="__main__":
    main()
