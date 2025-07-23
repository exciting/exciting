from os.path import dirname

import numpy as np

TUTORIAL_RUNDIR = "../run_tutorial_xc_functionals"


def test_xc_functionals():
    energy_vs_volume_GGA_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/GGA_PBE_SOL/energy-vs-volume")
    energy_vs_volume_GGA_reference = np.array([ [231.5687158849, -579.00840951],
                                              [238.9586580168, -579.01285013],
                                              [246.5041721916, -579.01616098],
                                              [254.2068789516, -579.01843820],
                                              [262.0683988388, -579.01976946],
                                              [270.0911626671, -579.02024024],
                                              [278.2743601632, -579.01992827],
                                              [286.6220426847, -579.01890781],
                                              [295.1350205018, -579.01724619],
                                              [303.8149141567, -579.01500868],
                                              [312.6633441916, -579.01225445]])

    assert np.allclose(energy_vs_volume_GGA_results, energy_vs_volume_GGA_reference,  atol=1.0e-8),(
                       "GGA energy vs volume values not equivalent to reference calculation")

    energy_vs_volume_LDA_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/LDA_PW/energy-vs-volume")
    energy_vs_volume_LDA_reference = np.array([[231.5687158849, -578.07161516],
                                              [238.9586580168, -578.07548055],
                                              [246.5041721916, -578.07824372],
                                              [254.2068789516, -578.08000041],
                                              [262.0683988388, -578.08084058],
                                              [270.0911626671, -578.08084512],
                                              [278.2743601632, -578.08009501],
                                              [286.6220426847, -578.07865992],
                                              [295.1350205018, -578.07661005],
                                              [303.8149141567, -578.07400722],
                                              [312.6633441916, -578.07091053]])

    assert np.allclose(energy_vs_volume_LDA_results, energy_vs_volume_LDA_reference, atol=1.0e-8),(
                       "LDA energy vs volume values not equivalent to reference calculation")

    energy_vs_volume_libxc_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/libxc_functional/energy-vs-volume")
    energy_vs_volume_libxc_reference = np.array([[231.5687158849, -580.05681318],
                                                [238.9586580168, -580.06184579],
                                                [246.5041721916, -580.06572902],
                                                [254.2068789516, -580.06856039],
                                                [262.0683988388, -580.07043173],
                                                [270.0911626671, -580.07142570],
                                                [278.2743601632, -580.07162449],
                                                [286.6220426847, -580.07109976],
                                                [295.1350205018, -580.06992224],
                                                [303.8149141567, -580.06815775],
                                                [312.6633441916, -580.06586460]])

    assert np.allclose(energy_vs_volume_libxc_results, energy_vs_volume_libxc_reference, atol=1.0e-8),(
                       "libxc energy vs volume values not equivalent to reference calculation")


def main():
    test_xc_functionals()

if __name__=="__main__":
    main()
