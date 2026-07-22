from os.path import dirname

import numpy as np

TUTORIAL_RUNDIR = "../run_diamond_elastic_strain/deformation-0"


def test_energy_vs_strain():
    energy_vs_strain_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/energy-vs-strain")
    energy_vs_strain_reference = np.array([[-0.10000000, -75.81585555],
                                            [-0.08000000, -75.84800729],
                                            [-0.06000000, -75.87014461],
                                            [-0.04000000, -75.88417959],
                                            [-0.02000000, -75.89175983],
                                            [ 0.00000000, -75.89408634],
                                            [ 0.02000000, -75.89221918],
                                            [ 0.04000000, -75.88701594],
                                            [ 0.06000000, -75.87910965],
                                            [ 0.08000000, -75.86908767],
                                            [ 0.10000000, -75.85741801]])

    assert np.allclose(energy_vs_strain_results, energy_vs_strain_reference, atol=1.0e-8), (
        "energy vs strain values not equivalent to reference calculation")


def main():
    test_energy_vs_strain()

if __name__=="__main__":
    main()
