"""
Tests for the gw_taskgroup_parser
"""

import numpy as np
import pytest

from excitingtools.exciting_dict_parsers.gw_taskgroup_parser import (
    parse_barc,
    parse_epsilon,
    parse_inverse_epsilon,
    parse_sgi,
)

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
    "matrix": np.array(
        [
            [complex(1.01e-4, -5.5e-8), complex(0.25, 0.77), complex(5.4e5, -1.1)],
            [complex(0.25, 0.77), complex(0.0, 0.0), complex(-5.4e5, 1.1)],
        ]
    )
}

square_matrix = """ 2
1 1 2 2
(1.01E-4,-5.5E-8) (0.000000000000000E+000,0.000000000000000E+000)
(5.4E+005,-1.1)
(-5.4E+005,1.1)
"""

reference_square_matrix = {
    "matrix": np.array([[complex(1.01e-4, -5.5e-8), complex(5.4e5, -1.1)], [complex(0.0, 0.0), complex(-5.4e5, 1.1)]])
}


@pytest.mark.parametrize(
    ["barc_file_str", "reference_barc"],
    [(rectangular_matrix, reference_rectangular_matrix), (square_matrix, reference_square_matrix)],
)
def test_parse_barc(barc_file_str, reference_barc, tmp_path):
    barc_file_path = tmp_path / "BARC_1.OUT"
    barc_file_path.write_text(barc_file_str)
    barc = parse_barc(barc_file_path.as_posix())
    A = reference_barc["matrix"]
    ref = {"CoulombMatrix": np.matmul(A.T.conj(), A)}
    np.testing.assert_allclose(barc["CoulombMatrix"], ref["CoulombMatrix"])


@pytest.mark.parametrize(
    ["sgi_file_str", "reference_sgi"],
    [(rectangular_matrix, reference_rectangular_matrix), (square_matrix, reference_square_matrix)],
)
def test_parse_sgi(sgi_file_str, reference_sgi, tmp_path):
    sgi_file_path = tmp_path / "SGI_1.OUT"
    sgi_file_path.write_text(sgi_file_str)
    sgi = parse_sgi(sgi_file_path.as_posix())
    A = reference_sgi["matrix"]
    ref = {"OverlapMatrix": np.matmul(A.T.conj(), A)}
    np.testing.assert_allclose(sgi["OverlapMatrix"], ref["OverlapMatrix"])


array_of_rank_3_example_1 = """ 3
1 1 1 2 2 1
(1.01E-4,-5.5E-8) (0.25,0.77)
(0.35,0.88) (0.000000000000000E+000,0.000000000000000E+000)
"""

reference_array_of_rank_3_example_1 = {
    "array": np.array(
        [[[complex(1.01e-4, -5.5e-8)], [complex(0.35, 0.88)]], [[complex(0.25, 0.77)], [complex(0.0, 0.0)]]]
    )
}

array_of_rank_3_example_2 = """ 3
1 1 1 3 4 2
(1.01E-4,-5.5E-8) (0.25,0.77)
(0.35,0.88) (0.000000000000000E+000,0.000000000000000E+000)
(-1.01E+004,-5.4E-8) (0.19,0.21) (1.5E-3,-9E6)
(0.07,0.00) (0.08,0.09) (1.1,2.2) (-4.5,-2.1) (-6.0E+000,8.0E-006)
(0.1,0.2) (1.5,2.5) (3.5,7.8) (7.0,8.7) (5.0,6.5) (8.0,9.1)
(0.01,0.00) (0.00,0.01) (0.02,0.00) (0.00,0.02) (0.3,0.0) (0.0,0.3)
"""

reference_array_of_rank_3_example_2 = {
    "array": np.array(
        [
            [
                [complex(1.01e-4, -5.5e-8), complex(0.1, 0.2)],
                [complex(0.0, 0.0), complex(7.0, 8.7)],
                [complex(1.5e-3, -9e6), complex(0.01, 0.00)],
                [complex(1.1, 2.2), complex(0.00, 0.02)],
            ],
            [
                [complex(0.25, 0.77), complex(1.5, 2.5)],
                [complex(-1.01e4, -5.4e-8), complex(5.0, 6.5)],
                [complex(0.07, 0.00), complex(0.00, 0.01)],
                [complex(-4.5, -2.1), complex(0.3, 0.0)],
            ],
            [
                [complex(0.35, 0.88), complex(3.5, 7.8)],
                [complex(0.19, 0.21), complex(8.0, 9.1)],
                [complex(0.08, 0.09), complex(0.02, 0.00)],
                [complex(-6, 8e-6), complex(0.0, 0.3)],
            ],
        ]
    )
}


@pytest.mark.parametrize(
    ["file_epsilon_str", "reference_epsilon"],
    [
        (array_of_rank_3_example_1, reference_array_of_rank_3_example_1),
        (array_of_rank_3_example_2, reference_array_of_rank_3_example_2),
    ],
)
def test_parse_epsilon(file_epsilon_str, reference_epsilon, tmp_path):
    epsilon_file_path = tmp_path / "EPSILON-GW_Q1.OUT"
    epsilon_file_path.write_text(file_epsilon_str)
    epsilon = parse_epsilon(epsilon_file_path.as_posix())
    np.testing.assert_allclose(epsilon["epsilon_tensor"], reference_epsilon["array"])


@pytest.mark.parametrize(
    ["file_inverse_epsilon_str", "reference_inverse_epsilon"],
    [
        (array_of_rank_3_example_1, reference_array_of_rank_3_example_1),
        (array_of_rank_3_example_2, reference_array_of_rank_3_example_2),
    ],
)
def test_parse_inverse_epsilon(file_inverse_epsilon_str, reference_inverse_epsilon, tmp_path):
    inverse_epsilon_file_path = tmp_path / "INVERSE-EPSILON_Q1.OUT"
    inverse_epsilon_file_path.write_text(file_inverse_epsilon_str)
    inverse_epsilon = parse_inverse_epsilon(inverse_epsilon_file_path.as_posix())
    np.testing.assert_allclose(inverse_epsilon["inverse_epsilon_tensor"], reference_inverse_epsilon["array"])
