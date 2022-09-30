""" Test for get_index function in eigenstates.py
"""

from excitingtools.eigenstates.eigenstates import get_k_point_index
import pytest
import numpy as np
import re


def test_get_k_point_index():
    k_points = np.array([[1.000000, 0.000000, 0.000000],
                         [0.988281, 0.011719, 0.000000],
                         [0.988281, 0.011719, 0.000000],
                         [0.976562, 0.023438, 0.000000],
                         [0.953125, 0.046875, 0.000000]])

    assert get_k_point_index([1.000000, 0.000000, 0.000000], k_points) == 0
    assert get_k_point_index([0.953125, 0.046875, 0.000000], k_points) == 4
    assert np.isnan(get_k_point_index([0.000000, 0.500000, 0.500000], k_points)), "No k-point matched"
    with pytest.raises(ValueError, match=re.escape("Found degenerate k-points at [[1], [2]]")):
        get_k_point_index([0.988281, 0.011719, 0.000000], k_points), "ValueError is returned for degenerate k-points"
