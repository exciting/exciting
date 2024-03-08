"""
Tests for the gw_taskgroup_parser
"""

import numpy as np
import pytest

from excitingtools.exciting_dict_parsers.gw_taskgroup_parser import parse_barc, \
    parse_sgi

rectangular_matrix = """ 2
1 1 2 3
(1.01E-4,-5.5E-8)
(0.25,0.77)
(0.25,0.77)
(0.000000000000000E+000,0.000000000000000E+000)
(5.4E+005,-1.1)
(-5.4E+005,1.1)
"""

reference_rectangular_matrix = {
    "matrix": np.array([
        [complex(1.01e-4, -5.5e-8), complex(0.25, 0.77), complex(5.4e5, -1.1)],
        [complex(0.25, 0.77), complex(0., 0.), complex(-5.4e5, 1.1)]
    ])
}

square_matrix = """ 2
1 1 2 2
(1.01E-4,-5.5E-8) (0.000000000000000E+000,0.000000000000000E+000)
(5.4E+005,-1.1)
(-5.4E+005,1.1)
"""

reference_square_matrix = {
    "matrix": np.array([
        [complex(1.01e-4, -5.5e-8), complex(5.4e5, -1.1)],
        [complex(0., 0.), complex(-5.4e5, 1.1)]
    ])
}


@pytest.mark.parametrize(["barc_file_str", "reference_barc"],
                         [(rectangular_matrix, reference_rectangular_matrix),
                          (square_matrix, reference_square_matrix)])
def test_parse_barc(barc_file_str, reference_barc, tmp_path):
    barc_file_path = tmp_path / "BARC_1.OUT"
    barc_file_path.write_text(barc_file_str)
    barc = parse_barc(barc_file_path.as_posix())
    A = reference_barc["matrix"]
    ref = {"CoulombMatrix": np.matmul(A.T.conj(), A)}
    assert np.allclose(barc["CoulombMatrix"], ref["CoulombMatrix"])


@pytest.mark.parametrize(["sgi_file_str", "reference_sgi"],
                         [(rectangular_matrix, reference_rectangular_matrix),
                          (square_matrix, reference_square_matrix)])
def test_parse_sgi(sgi_file_str, reference_sgi, tmp_path):
    sgi_file_path = tmp_path / "SGI_1.OUT"
    sgi_file_path.write_text(sgi_file_str)
    sgi = parse_sgi(sgi_file_path.as_posix())
    A = reference_sgi["matrix"]
    ref = {"OverlapMatrix": np.matmul(A.T.conj(), A)}
    assert np.allclose(sgi["OverlapMatrix"], ref["OverlapMatrix"])
