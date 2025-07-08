from os.path import dirname

import numpy as np
from excitingtools.exciting_dict_parsers.bse_parser import parse_EPSILON_NAR

TUTORIAL_RUNDIR = "../run_LiF_TDDFT_kernels"
REFERENCE_DIR = "reference_many_body_kernels_for_tddft_calculations"


def test_many_body_kernels_tddft(file_name: str):
    """Automatically test results of many_body_kernels_for_tddft_calculations tutorial.
    """

    epsilon_reference = parse_EPSILON_NAR(f"{dirname(__file__)}/{REFERENCE_DIR}/{file_name}")
    epsilon_result = parse_EPSILON_NAR(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/{file_name}")

    assert np.allclose(epsilon_result['frequency'],
                       epsilon_reference['frequency'], atol = 1e-6), (
        "Freqeuncy grid not equivalent \
         to reference calculation")

    assert np.allclose(epsilon_result['real_oscillator_strength'],
                       epsilon_reference['real_oscillator_strength'], atol = 1e-3), (
        "Real part of dielectric function not equivalent \
        to reference calculation")

    assert np.allclose(epsilon_result['imag_oscillator_strength'],
                       epsilon_reference['imag_oscillator_strength'], atol = 1e-4), (
        "Imaginary part of dielectric function not equivalent \
        to reference calculation")

    assert np.allclose(epsilon_result['real_oscillator_strength_kkt'],
                       epsilon_reference['real_oscillator_strength_kkt'], atol = 1e-3), (
        "Real part of dielectric function (by Kramers-Kronig) \
         not equivalent to reference calculation")


def main():
    epsilon_files = ["BSE_kernel/EPSILON_NAR_FXCMB1_OC11_QMT001.OUT", "LRC_RBO_kernel/EPSILON_FXCRBO_OC11_QMT001.OUT",
                     "LRC_RBO_kernel/EPSILON_NLF_FXCLRCstatic_OC11_QMT001.OUT"]

    for file in epsilon_files:
        test_many_body_kernels_tddft(file)


if __name__ == "__main__":
    main()
