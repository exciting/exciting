from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from excitingtools.exciting_dict_parsers.phonon_parser import (
    parse_dyn_out,
    parse_epsinf_out,
    parse_phonon_out,
    parse_zstar_out,
)
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


dyn_out_content = """  0.2515698486       0.000000000     : is =    1, ia =    1, ip =    1
   0.000000000       0.000000000     : is =    1, ia =    1, ip =    2
   0.000000000       0.000000000     : is =    1, ia =    1, ip =    3
 -0.2515698486       0.000000000     : is =    2, ia =    1, ip =    1
   0.000000000       0.000000000     : is =    2, ia =    1, ip =    2
   0.000000000       0.000000000     : is =    2, ia =    1, ip =    3
"""


def test_parse_dyn_out(tmp_path: Path) -> None:
    dyn_out_file = tmp_path / "DYN_Q0000_0000_0000_S01_A001_P1.OUT"
    dyn_out_file.write_text(dyn_out_content)

    dyn_out = parse_dyn_out(dyn_out_file)

    ref_dyn_out = {
        "1": {"species": 1, "atom": 1, "polarisation": 1, "dynmat_real": 0.2515698486, "dynmat_imag": 0.0},
        "2": {"species": 1, "atom": 1, "polarisation": 2, "dynmat_real": 0.0, "dynmat_imag": 0.0},
        "3": {"species": 1, "atom": 1, "polarisation": 3, "dynmat_real": 0.0, "dynmat_imag": 0.0},
        "4": {"species": 2, "atom": 1, "polarisation": 1, "dynmat_real": -0.2515698486, "dynmat_imag": 0.0},
        "5": {"species": 2, "atom": 1, "polarisation": 2, "dynmat_real": 0.0, "dynmat_imag": 0.0},
        "6": {"species": 2, "atom": 1, "polarisation": 3, "dynmat_real": 0.0, "dynmat_imag": 0.0},
    }

    assert len(dyn_out) == 6
    for key, value in ref_dyn_out.items():
        assert value == pytest.approx(dyn_out[key])


epsinf_content = """# High frequency dielectric tensor (clamped nuclei).
#
        4.5238287243        0.0000000000        0.0000000000
        0.0000000000        4.5238287243        0.0000000000
        0.0000000000        0.0000000000        4.5238287243
"""


def test_parse_epsinf(tmp_path: Path) -> None:
    epsinf_file = tmp_path / "EPSINF.OUT"
    epsinf_file.write_text(epsinf_content)

    epsinf_out = parse_epsinf_out(epsinf_file)
    assert len(epsinf_out) == 1
    epsinf_array = epsinf_out["epsinf"]

    epsinf_ref = np.array([[4.52382872, 0.0, 0.0], [0.0, 4.52382872, 0.0], [0.0, 0.0, 4.52382872]])
    assert_allclose(epsinf_array, epsinf_ref)


zstar_content = """# Born effective charge tensors for all atoms.
# Rows correspond to E-field direction.
# Columns correspond to atom displacement direction.
# Acoustic sum rule has been imposed.
#
# species  1 atom   1 (B  1) :      0.000000     0.000000     0.000000
        1.9552213065        0.0000000000        0.0000000000
        0.0000000000        1.9552213065        0.0000000000
        0.0000000000        0.0000000000        1.9552213065
# species  2 atom   1 (N  1) :      0.250000     0.250000     0.250000
       -1.9552213065        0.0000000000        0.0000000000
        0.0000000000       -1.9552213065        0.0000000000
        0.0000000000        0.0000000000       -1.9552213065
# Acoustic sum rule correction (add to each tensor above to get original value)
        0.0004145685        0.0000000000        0.0000000000
        0.0000000000        0.0004145685        0.0000000000
        0.0000000000        0.0000000000        0.0004145685
"""


def test_parse_zstar(tmp_path: Path) -> None:
    zstar_file = tmp_path / "ZSTAR.OUT"
    zstar_file.write_text(zstar_content)

    zstar_out = parse_zstar_out(zstar_file)

    assert len(zstar_out) == 2

    ref_zstar_out = {
        "atoms": {
            "1": {
                "species": 1,
                "atom": 1,
                "symbol": "B",
                "species_atom": 1,
                "position": np.array([0.0, 0.0, 0.0]),
                "tensor": np.array([[1.95522131, 0.0, 0.0], [0.0, 1.95522131, 0.0], [0.0, 0.0, 1.95522131]]),
            },
            "2": {
                "species": 2,
                "atom": 1,
                "symbol": "N",
                "species_atom": 1,
                "position": np.array([0.25, 0.25, 0.25]),
                "tensor": np.array([[-1.95522131, 0.0, 0.0], [0.0, -1.95522131, 0.0], [0.0, 0.0, -1.95522131]]),
            },
        },
        "acoustic_sum_rule_correction": np.array(
            [[0.0004145685, 0.0, 0.0], [0.0, 0.0004145685, 0.0], [0.0, 0.0, 0.0004145685]]
        ),
    }

    assert_allclose(zstar_out["acoustic_sum_rule_correction"], ref_zstar_out["acoustic_sum_rule_correction"])
    ref_atoms = ref_zstar_out["atoms"]
    atoms = zstar_out["atoms"]

    for key, value in ref_atoms.items():
        for key2, value2 in value.items():
            if isinstance(value2, np.ndarray):
                assert_allclose(value2, atoms[key][key2])
            else:
                assert value2 == atoms[key][key2]
