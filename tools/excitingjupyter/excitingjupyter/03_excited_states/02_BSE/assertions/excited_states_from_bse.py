from os.path import dirname

import numpy as np
from excitingtools import parse

TUTORIAL_EXCITED_STATES_FROM_BSE_RUNDIR = "../run_LiF_BSE/BSE/EPSILON"
REFERENCE_DIR = "reference_excited_states_from_BSE/EPSILON"

def test_tutorial_excited_states_from_bse(file_name: str):
    """Automatically test results of excited_states_from_bse tutorial, LiF bulk calculation.
    """
    
    epsilon_reference = parse(f"{dirname(__file__)}/{REFERENCE_DIR}/{file_name}")
    epsilon_result = parse(f"{dirname(__file__)}/{TUTORIAL_EXCITED_STATES_FROM_BSE_RUNDIR}/{file_name}")

    assert np.allclose(epsilon_result['frequency'],
                       epsilon_reference['frequency'], atol=2e-6), (
                       "Freqeuncy grid not equivalent \
                        to reference calculation")

    assert np.allclose(epsilon_result['real_oscillator_strength'],
                       epsilon_reference['real_oscillator_strength'], atol=2e-6), (
                        "Real part of dielectric function not equivalent \
                        to reference calculation")

    assert np.allclose(epsilon_result['imag_oscillator_strength'],
                       epsilon_reference['imag_oscillator_strength'], atol=2e-6), (
                        "Imaginary part of dielectric function not equivalent \
                        to reference calculation")
                        
    assert np.allclose(epsilon_result['real_oscillator_strength_kkt'],
                       epsilon_reference['real_oscillator_strength_kkt'], atol=2e-6), (
                        "Real part of dielectric function (by Kramers-Kronig) \
                         not equivalent to reference calculation")

def main():
    epsilon_files = ["EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT", "EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT",
                     "EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT"]

    for file in epsilon_files:
        test_tutorial_excited_states_from_bse(file)

if __name__=="__main__":
    main()
