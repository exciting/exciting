from os.path import dirname

import numpy as np
from excitingtools.exciting_dict_parsers.properties_parser import parse_band_structure_xml

TUTORIAL_RUNDIR = "../run_tutorial_exact_exchange_calculations"
REFERENCE_DIR = "reference_exact_exchange_calculations"

def test_bandstructure_exact_exchange():
    bandstructure_LDA_reference = parse_band_structure_xml(
        f"{dirname(__file__)}/{REFERENCE_DIR}/bandstructure_LDA.xml")
    bandstructure_LDA_results = parse_band_structure_xml(
        f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/LDA/bandstructure.xml")

    assert np.allclose(bandstructure_LDA_results["k_points_along_band"],
                       bandstructure_LDA_reference["k_points_along_band"], atol=1.0e-6),(
                       "LDA k-point values not equivalent to reference calculation")
    assert np.allclose(bandstructure_LDA_results["band_energies"],
                       bandstructure_LDA_reference["band_energies"], atol=1.0e-6),(
                       "LDA band energies not equivalent to reference calculation")

    bandstructure_EXX_reference = parse_band_structure_xml(
        f"{dirname(__file__)}/{REFERENCE_DIR}/bandstructure_EXX.xml")
    bandstructure_EXX_results = parse_band_structure_xml(
        f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/EXX/bandstructure.xml")

    assert np.allclose(bandstructure_EXX_results["k_points_along_band"],
                       bandstructure_EXX_reference["k_points_along_band"], atol=1.0e-6),(
                       "EXX k-point values not equivalent to reference calculation")
    assert np.allclose(bandstructure_EXX_results["band_energies"],
                       bandstructure_EXX_reference["band_energies"], atol=1.0e-6),(
                       "EXX band energies not equivalent to reference calculation")


def main():
    test_bandstructure_exact_exchange()

if __name__=="__main__":
    main()
