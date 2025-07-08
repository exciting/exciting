from os.path import dirname

import numpy as np
from excitingtools.exciting_dict_parsers.bse_parser import parse_EPSILON_NAR, parse_LOSS_NAR

TUTORIAL_RUNDIR = "../run_Ag_q_tddft/q001_test_1"
REFERENCE_DIR = "reference_q_dependent_tddft"


def test_dielectric_q_tddft(file_name: str):
    """Test results of dielectric function calculations for fastBSE.
    """

    epsilon_reference = parse_EPSILON_NAR(f"{dirname(__file__)}/{REFERENCE_DIR}/{file_name}")
    epsilon_results = parse_EPSILON_NAR(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/{file_name}")

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

def test_loss_q_tddft(file_name: str):
    """Test results of loss function calculations for fastBSE.
    """

    loss_reference = parse_LOSS_NAR(f"{dirname(__file__)}/{REFERENCE_DIR}/{file_name}")
    loss_results = parse_LOSS_NAR(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/{file_name}")

    assert np.allclose(loss_results['frequency'], loss_reference['frequency']), \
        f"Frequency grid not equivalent to reference calculation."

    assert np.allclose(loss_results['real_oscillator_strength'], loss_reference['real_oscillator_strength']), \
        f"Real part of loss function not equivalent to reference calculation."

    assert np.allclose(loss_results['imag_oscillator_strength'], loss_reference['imag_oscillator_strength']), \
        f"Imaginary part of loss function not equivalent to reference calculation."

def main():
    test_dielectric_q_tddft("EPSILON_FXCRPA_OC11_QMT001.OUT")

    test_loss_q_tddft("LOSS_FXCRPA_OC11_QMT001.OUT")

if __name__=="__main__":
    main()
