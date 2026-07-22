from os.path import dirname

import numpy as np
from excitingtools.exciting_dict_parsers.properties_parser import parse_band_structure_xml

TUTORIAL_RUNDIR = "../run_spin_tutorial"
REFERENCE_DIR = "reference_spin_polarized_calculations"

def test_bandstructure_spin_polarized():
    bandstructure_spin_reference = parse_band_structure_xml(
        f"{dirname(__file__)}/{REFERENCE_DIR}/bandstructure_spin.xml")
    bandstructure_spin_results = parse_band_structure_xml(
        f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/SPIN/bandstructure.xml")

    assert np.allclose(bandstructure_spin_results["k_points_along_band"],
                       bandstructure_spin_reference["k_points_along_band"], atol=1.0e-6),(
                       "Spin k-point values not equivalent to reference calculation")
    assert np.allclose(bandstructure_spin_results["band_energies"],
                       bandstructure_spin_reference["band_energies"], atol=1.0e-6),(
                       "Spin band energies not equivalent to reference calculation")

    bandstructure_afm_reference = parse_band_structure_xml(
        f"{dirname(__file__)}/{REFERENCE_DIR}/bandstructure_afm.xml")
    bandstructure_afm_results = parse_band_structure_xml(
        f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/AFM-B_4.0/bandstructure.xml")

    assert np.allclose(bandstructure_afm_results["k_points_along_band"],
                       bandstructure_afm_reference["k_points_along_band"], atol=1.0e-6),(
                       "AFM k-point values not equivalent to reference calculation")
    assert np.allclose(bandstructure_afm_results["band_energies"],
                       bandstructure_afm_reference["band_energies"], atol=1.0e-6),(
                       "AFM band energies not equivalent to reference calculation")

def test_dos_spin_polarized():
    dos_spin_reference = np.loadtxt(f"{dirname(__file__)}/{REFERENCE_DIR}/TDOS_spin.OUT")
    dos_spin_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/SPIN/TDOS.OUT")

    rel_error = np.linalg.norm(dos_spin_results-dos_spin_reference)/np.linalg.norm(dos_spin_reference)
    assert rel_error <= 1e-5,(
                       f"Spin DOS not equivalent to reference calculation. Max err: {abs(dos_spin_results - dos_spin_reference).max()}, Rel err: {rel_error}")

    dos_afm_reference = np.loadtxt(f"{dirname(__file__)}/{REFERENCE_DIR}/TDOS_afm.OUT")
    dos_afm_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/AFM-B_4.0/TDOS.OUT")

    assert np.allclose(dos_afm_results, dos_afm_reference, atol=1.0e-6),(
                       "AFM DOS not equivalent to reference calculation")

def main():
    test_bandstructure_spin_polarized()
    test_dos_spin_polarized()

if __name__=="__main__":
    main()
