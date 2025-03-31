from os.path import dirname

import numpy as np

TUTORIAL_RUNDIR = "../run_diamond_phonons"
REFERENCE_DIR = "reference_lattice_dynamics_of_diamond_and_zincblende_structure_crystals"

def test_phonon_dos(file_name: str):
    dos_reference = np.loadtxt(f"{dirname(__file__)}/{REFERENCE_DIR}/{file_name}")
    dos_results = np.loadtxt(f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/{file_name}")

    assert np.allclose(dos_results, dos_reference, atol=1.0e-8),(
                       "Phonon DOS not equivalent to reference calculation")

def main():
    test_phonon_dos("PHDOS.OUT")

if __name__=="__main__":
    main()
