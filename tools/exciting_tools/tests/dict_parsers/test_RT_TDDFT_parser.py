"""
Test for the RT_TDDFT_parser
"""

import numpy as np
import pytest

from excitingtools.exciting_dict_parsers.RT_TDDFT_parser import parse_occupations, parse_proj_screenshots

occupations_str = """ik =       1
      1     1.97843104
      2     1.10327541
      3     0.57293646
      4     0.31164780
      5     0.21001308
      6     0.00093199

ik =       2
      1     1.99998297
      2     1.99990887
      3     0.00131209
      4     0.00146791
      5     0.00079852
      6     0.00095539

"""

occ_ref = {
    "ik": [1, 2],
    "occupations": [
        np.array([1.97843104, 1.10327541, 0.57293646, 0.31164780, 0.21001308, 0.00093199]),
        np.array([1.99998297, 1.99990887, 0.00131209, 0.00146791, 0.00079852, 0.00095539]),
    ],
}


def test_parse_occupations(tmp_path) -> None:
    file_path = tmp_path / "OCCSV_TXT_10.OUT"
    file_path.write_text(occupations_str)
    occ = parse_occupations(file_path.as_posix())
    assert occ["ik"] == occ_ref["ik"]
    np.testing.assert_allclose(occ["occupations"], occ_ref["occupations"])


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
def test_parse_proj_screenshots(proj_file_str, reference_parsed_dict, tmp_path) -> None:
    proj_file_path = tmp_path / "PROJ_0.OUT"
    proj_file_path.write_text(proj_file_str)
    proj_out = parse_proj_screenshots(proj_file_path.as_posix())
    is_equal = proj_out["ik"] == reference_parsed_dict["ik"]
    key = "projection"
    is_equal = is_equal and all([np.allclose(x, y) for (x, y) in zip(proj_out[key], reference_parsed_dict[key])])
    assert is_equal
