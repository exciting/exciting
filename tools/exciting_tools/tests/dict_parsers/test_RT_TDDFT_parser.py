"""
Test for the RT_TDDFT_parser
"""

import numpy as np
import pytest

from excitingtools.exciting_dict_parsers.RT_TDDFT_parser import parse_proj_screenshots

proj_file_str_square_matrices = """ ik:          1
   1.00000   0.00000
   0.00000   1.00000
 ik:          2
   1.00000   0.00000
   0.00000   1.00000
"""

reference_parsed_proj_square_matrices = {
    "ik": [1, 2],
    "projection": [np.array([[1.0, 0.0], [0.0, 1.0]]), np.array([[1.0, 0.0], [0.0, 1.0]])],
}

proj_file_str_rectangular_matrices = """ ik:          1
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000
 ik:          2
   0.60000   0.80000   0.00000
   0.00000   0.00000   1.00000
"""

reference_parsed_proj_rectangular_matrices = {
    "ik": [1, 2],
    "projection": [np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]), np.array([[0.6, 0.8, 0.0], [0.0, 0.0, 1.0]])],
}


@pytest.mark.parametrize(
    ["proj_file_str", "reference_parsed_dict"],
    [
        (proj_file_str_square_matrices, reference_parsed_proj_square_matrices),
        (proj_file_str_rectangular_matrices, reference_parsed_proj_rectangular_matrices),
    ],
)
def test_parse_proj_screenshots(proj_file_str, reference_parsed_dict, tmp_path):
    proj_file_path = tmp_path / "PROJ_0.OUT"
    proj_file_path.write_text(proj_file_str)
    proj_out = parse_proj_screenshots(proj_file_path.as_posix())
    is_equal = proj_out["ik"] == reference_parsed_dict["ik"]
    key = "projection"
    is_equal = is_equal and all([np.allclose(x, y) for (x, y) in zip(proj_out[key], reference_parsed_dict[key])])
    assert is_equal
