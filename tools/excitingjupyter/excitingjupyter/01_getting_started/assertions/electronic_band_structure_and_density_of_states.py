from os.path import dirname

import numpy as np
from excitingtools.exciting_dict_parsers.properties_parser import parse_band_structure_xml

TUTORIAL_RUNDIR = "../run_bs_dos"
REFERENCE_DIR = "reference_electronic_band_structure_and_density_of_states"

def test_bandstructure(file_name: str):
    bandstructure_reference = parse_band_structure_xml(f"{dirname(__file__)}/{REFERENCE_DIR}/{file_name}")
    bandstructure_results = parse_band_structure_xml(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/{file_name}")

    assert np.allclose(bandstructure_results["k_points_along_band"],
                       bandstructure_reference["k_points_along_band"], atol=1.0e-8),(
                       "k-point values not equivalent to reference calculation")
    assert np.allclose(bandstructure_results["band_energies"],
                       bandstructure_reference["band_energies"], atol=1.0e-8),(
                       "band energies not equivalent to reference calculation")

def test_dos(file_name: str):
    dos_reference = np.loadtxt(f"{dirname(__file__)}/{REFERENCE_DIR}/{file_name}")
    dos_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/{file_name}")

    assert np.allclose(dos_results, dos_reference, atol=1.0e-8),(
                       "DOS not equivalent to reference calculation")

def main():
    test_bandstructure("bandstructure.xml")
    test_dos("TDOS.OUT")

if __name__=="__main__":
    main()
