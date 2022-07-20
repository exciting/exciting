import numpy as np
import os.path
from excitingtools import parser_chooser


def test_tutorial_excited_states_from_bse(epsilon_result: dict, path_to_reference: str):
    """Automatically test results of tutorial_excited_states_from_bse notebook, LiF bulk calculation.
    """
    
    epsilon_reference = parser_chooser(os.path.join(path_to_reference,"EPSILON/EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT"))

    assert np.allclose(epsilon_result['frequency'],
                       epsilon_reference['frequency']), (
                       "Freqeuncy grid not equivalent \
                        to reference calculation")

    assert np.allclose(epsilon_result['real_oscillator_strength'],
                       epsilon_reference['real_oscillator_strength']), (
                        "Real part of dielectric function not equivalent \
                        to reference calculation")

    assert np.allclose(epsilon_result['imag_oscillator_strength'],
                       epsilon_reference['imag_oscillator_strength']), (
                        "Imaginary part of dielectric function not equivalent \
                        to reference calculation")
                        
    assert np.allclose(epsilon_result['real_oscillator_strength_kkt'],
                       epsilon_reference['real_oscillator_strength_kkt']), (
                        "Real part of dielectric function (by Kramers-Kronig) \
                         not equivalent to reference calculation")