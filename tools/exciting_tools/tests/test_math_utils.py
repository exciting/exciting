"""
Tests for functions in math_utils.py
"""
import numpy as np

from excitingtools.math_utils import unit_vector


def test_unit_vector():
    x = np.array([0.23, 0.32, 1.394, 99.])
    norm_x = np.linalg.norm(unit_vector(x))
    assert np.isclose(norm_x, 1., atol=1.e-8)

