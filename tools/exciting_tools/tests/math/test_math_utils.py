"""
Tests for functions in math_utils.py
"""
import numpy as np

from excitingtools.math.math_utils import unit_vector


def test_unit_vector():
    """Test unit_vector returns a vector with norm of one."""
    x = np.array([0.23, 0.32, 1.394, 99.])
    norm_x = np.linalg.norm(unit_vector(x))
    assert np.isclose(norm_x, 1., atol=1.e-8)
