from os.path import dirname

from excitingtools.exciting_dict_parsers.phonon_parser import parse_phonon_out

TUTORIAL_RUNDIR = "../run_diamond_phonons"
REFERENCE_DIR = "reference_lattice_dynamics_of_diamond_and_zincblende_structure_crystals"

def test_phonon_data():
    phonon_reference = parse_phonon_out(f"{dirname(__file__)}/{REFERENCE_DIR}/PHONON.OUT")
    phonon_results = parse_phonon_out( f"{dirname(__file__)}/{TUTORIAL_RUNDIR}/PHONON.OUT")

    assert phonon_results == phonon_reference,(
                       "PHONON.OUT data not equivalent to reference calculation")


def main():
    test_phonon_data()

if __name__=="__main__":
    main()
