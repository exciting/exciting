from pathlib import Path

import pytest

from excitingtools.exciting_dict_parsers.phonon_parser import parse_phonon_out
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def phonon_out_mock(tmp_path: Path) -> MockFile:
    """Mock PHONON.OUT data with two q-points and corresponding modes."""
    phonon_out_str = """
     1   0.000000000       0.000000000       0.000000000     : q-point, vpl

     1  0.8231806350E-10 : mode, frequency
   1   1   1  0.7071067691       0.000000000     : species, atom, polarisation, eigenvector
   1   1   2   0.000000000       0.000000000    
   1   1   3   0.000000000       0.000000000    
   1   2   1  0.7071067691       0.000000000    
   1   2   2   0.000000000       0.000000000    
   1   2   3   0.000000000       0.000000000    

     2  0.8231806350E-10 : mode, frequency
   1   1   1   0.000000000       0.000000000     : species, atom, polarisation, eigenvector
   1   1   2  0.7071067691       0.000000000    
   1   1   3   0.000000000       0.000000000    
   1   2   1   0.000000000       0.000000000    
   1   2   2  0.7071067691       0.000000000    
   1   2   3   0.000000000       0.000000000    

     2  0.5000000000      0.5000000000       0.000000000     : q-point, vpl

     1  0.3567275615E-02 : mode, frequency
   1   1   1   0.000000000       0.000000000     : species, atom, polarisation, eigenvector
   1   1   2   0.000000000       0.000000000    
   1   1   3  0.7071067691       0.000000000    
   1   2   1   0.000000000       0.000000000    
   1   2   2  0.7071067691       0.000000000    
   1   2   3   0.000000000       0.000000000    

     2  0.3567275615E-02 : mode, frequency
   1   1   1   0.000000000       0.000000000     : species, atom, polarisation, eigenvector
   1   1   2  0.7071067691       0.000000000    
   1   1   3   0.000000000       0.000000000    
   1   2   1   0.000000000       0.000000000    
   1   2   2   0.000000000       0.000000000    
   1   2   3  0.7071067691       0.000000000    
    """
    phonon_out_file = tmp_path / "PHONON.OUT"
    phonon_out_file.write_text(phonon_out_str)
    return MockFile(phonon_out_file, phonon_out_str)


def test_parse_phonon_out(phonon_out_mock: MockFile) -> None:
    phonon_data = parse_phonon_out(phonon_out_mock.file)

    assert len(phonon_data) == 2

    q1 = phonon_data["1"]
    assert q1["q_vector"] == [0.000000000, 0.000000000, 0.000000000]
    assert len(q1["modes"]) == 2

    mode1_q1 = q1["modes"][0]
    assert mode1_q1["mode_index"] == "1"
    assert mode1_q1["frequency"] == 0.8231806350e-10
    assert len(mode1_q1["eigenvector_info"]) == 6

    q2 = phonon_data["2"]
    assert q2["q_vector"] == [0.5000000000, 0.5000000000, 0.000000000]
    assert len(q2["modes"]) == 2

    mode1_q2 = q2["modes"][0]
    assert mode1_q2["mode_index"] == "1"
    assert mode1_q2["frequency"] == 0.3567275615e-02
    assert len(mode1_q2["eigenvector_info"]) == 6

    eig1 = mode1_q1["eigenvector_info"][0]
    assert eig1["species"] == 1
    assert eig1["atom"] == 1
    assert eig1["polarisation"] == 1
    assert eig1["eigenvector_component_real"] == 0.7071067691
    assert eig1["eigenvector_component_imag"] == 0.000000000
