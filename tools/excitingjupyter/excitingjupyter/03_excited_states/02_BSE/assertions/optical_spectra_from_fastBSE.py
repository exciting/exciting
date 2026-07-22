from os.path import dirname

import numpy as np
from excitingtools.exciting_dict_parsers.bse_parser import parse_EPSILON_NAR, parse_LOSS_NAR

TUTORIAL_RUNDIR = "../run_fastBSE_diamond"
REFERENCE_DIR = "reference_optical_spectra_from_fastBSE"


def test_dielectric_fastBSE(file_name: str):
    """Test results of dielectric function calculations for fastBSE.
    """

    epsilon_reference = parse_EPSILON_NAR(f"{dirname(__file__)}/{REFERENCE_DIR}/{file_name}")
    epsilon_results = parse_EPSILON_NAR(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/EPSILON/{file_name}")

    assert np.allclose(epsilon_results['frequency'],
                       epsilon_reference['frequency']), (
        "Freqeuncy grid not equivalent \
         to reference calculation")

    assert np.allclose(epsilon_results['real_oscillator_strength'],
                       epsilon_reference['real_oscillator_strength']), (
        "Real part of dielectric function not equivalent \
        to reference calculation")

    assert np.allclose(epsilon_results['imag_oscillator_strength'],
                       epsilon_reference['imag_oscillator_strength']), (
        "Imaginary part of dielectric function not equivalent \
        to reference calculation")

    assert np.allclose(epsilon_results['real_oscillator_strength_kkt'],
                       epsilon_reference['real_oscillator_strength_kkt']), (
        "Real part of dielectric function (by Kramers-Kronig) \
         not equivalent to reference calculation")

def test_loss_fastBSE(file_name: str):
    """Test results of loss function calculations for fastBSE.
    """

    loss_reference = parse_LOSS_NAR(f"{dirname(__file__)}/{REFERENCE_DIR}/{file_name}")
    loss_results = parse_LOSS_NAR(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/LOSS/{file_name}")

    assert np.allclose(loss_results['frequency'], loss_reference['frequency']), \
        f"Frequency grid not equivalent to reference calculation."

    assert np.allclose(loss_results['real_oscillator_strength'], loss_reference['real_oscillator_strength']), \
        f"Real part of loss function not equivalent to reference calculation."

    assert np.allclose(loss_results['imag_oscillator_strength'], loss_reference['imag_oscillator_strength']), \
        f"Imaginary part of loss function not equivalent to reference calculation."

def main():
    test_dielectric_fastBSE("EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT")

    test_loss_fastBSE("LOSS_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT")

if __name__=="__main__":
    main()
