from os.path import dirname

import numpy as np
from excitingtools.exciting_dict_parsers.bse_parser import parse_EXCITON_NAR_BSE


TUTORIAL_RUNDIR = "../run_LiF_exciton/EXCITON"
REFERENCE_DIR = "reference_exciton_analysis_and_visualization"


def test_exciton(file_name: str):
    """Automatically test results of tutorial_excited_states_from_bse notebook, LiF bulk calculation.
    """

    exciton_reference = parse_EXCITON_NAR_BSE(f"{dirname(__file__)}/{REFERENCE_DIR}/{file_name}")
    exciton_result = parse_EXCITON_NAR_BSE(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/{file_name}")

    assert np.allclose(exciton_result['energy'],
                       exciton_reference['energy']), (
        "Excitation energies not equivalent \
         to reference calculation")

    assert np.allclose(exciton_result['energy_shifted'],
                       exciton_reference['energy_shifted'], atol=1e-5), (
        "Exciton binding energies not equivalent \
        to reference calculation")

    assert np.allclose(exciton_result['abs_oscillator_strength'],
                       exciton_reference['abs_oscillator_strength'], atol=1e-5), (
        "Oscillator strengths not equivalent \
        to reference calculation")

    assert np.allclose(exciton_result['real_oscillator_strength'],
                       exciton_reference['real_oscillator_strength'], atol=1e-5), (
        "Real part of oscillator strengths \
         not equivalent to reference calculation")

    assert np.allclose(exciton_result['imaginary_oscillator_strength'],
                       exciton_reference['imaginary_oscillator_strength'], atol=1e-5), (
        "Imaginary part of oscillator strengths \
         not equivalent to reference calculation")

def main():
    test_exciton("EXCITON_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT")

if __name__ == "__main__":
    main()