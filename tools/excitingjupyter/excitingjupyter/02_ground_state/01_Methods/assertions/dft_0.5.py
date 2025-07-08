from os.path import dirname

import numpy as np

TUTORIAL_RUNDIR = "../run_tutorial_LDA_0.5/cut-Si"


def test_dft_half():
    bandgap_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/bandgap-with-all")
    bandgap_reference = np.array([[0.0, 0.6059858708],
                    [0.5, 0.6191368617],
                    [1.0, 0.6226174701],
                    [1.5, 0.6392923355],
                    [2.0, 0.7316698216],
                    [2.5, 0.9076733399],
                    [3.0, 1.1026573444],
                    [3.5, 1.1982754346],
                    [3.6, 1.2020510144],
                    [3.7, 1.2008270463],
                    [3.8, 1.1947749618],
                    [3.9, 1.1841189830],
                    [4.0, 1.1691102708],
                    [4.5, 1.0394795805],
                    [5.0, 0.8522698691]])

    assert np.allclose(bandgap_results, bandgap_reference,  atol=1.0e-7), (
                    "bandgap values not equivalent to reference calculation")


def main():
    test_dft_half()

if __name__=="__main__":
    main()
