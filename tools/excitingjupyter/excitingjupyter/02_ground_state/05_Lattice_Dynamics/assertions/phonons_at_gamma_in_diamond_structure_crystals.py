from os.path import dirname

import numpy as np

TUTORIAL_RUNDIR = "../run_diamond_phonon_g/gamma_1"


def test_phonons_gamma():
    energy_vs_displacement_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/energy-vs-displacement")
    energy_vs_displacement_reference = np.array([[-0.02500000, -75.9459459514],
                                                [-0.02000000, -75.9570243210],
                                                [-0.01500000, -75.9651916320],
                                                [-0.01000000, -75.9707282135],
                                                [-0.00500000, -75.9738857721],
                                                [ 0.00000100, -75.9748914397],
                                                [ 0.00500000, -75.9739487117],
                                                [ 0.01000000, -75.9712392120],
                                                [ 0.01500000, -75.9669252222],
                                                [ 0.02000000, -75.9611524361],
                                                [ 0.02500000, -75.9540536437]])

    assert np.allclose(energy_vs_displacement_results, energy_vs_displacement_reference,  atol=1.0e-8),(
                       "energy vs displacement values not equivalent to reference calculation")

    energy_vs_force_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/energy-vs-force")
    energy_vs_force_reference = np.array([[-0.02500000,  0.1244712479],
                                        [-0.02000000,  0.0943003992],
                                        [-0.01500000,  0.0670517692],
                                        [-0.01000000,  0.0424301316],
                                        [-0.00500000,  0.0201701688],
                                        [ 0.00000100,  0.0000278896],
                                        [ 0.00500000, -0.0182042429],
                                        [ 0.01000000, -0.0347361620],
                                        [ 0.01500000, -0.0497406137],
                                        [ 0.02000000, -0.0633704287],
                                        [ 0.02500000, -0.0757527934]])

    assert np.allclose(energy_vs_force_results, energy_vs_force_reference, atol=1.0e-8),(
                       "energy vs force values not equivalent to reference calculation")



def main():
    test_phonons_gamma()

if __name__=="__main__":
    main()
