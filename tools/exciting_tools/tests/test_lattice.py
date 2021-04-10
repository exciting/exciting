"""
Tests for functions in lattice.py 
"""
import numpy as np

from ..exciting_tools.lattice import reciprocal_lattice_vectors, parallelpiped_volume
from ..exciting_tools.math_utils import triple_product


def test_parallelpiped_volume():

    # FCC lattice vectors with arbitrary lattice constant
    a = 3.2 * np.array([[0., 1., 1.],
                        [1., 0., 1.],
                        [1., 1., 0.]])

    triple_product_result = triple_product(a[:, 0], a[:, 1], a[:, 2])
    volume = parallelpiped_volume(np.array([a[:, 0], a[:, 1], a[:, 2]]))

    assert np.isclose(triple_product_result, 65.5360)
    assert np.isclose(volume, 65.5360)

    triple_product_result = triple_product(a[:, 0], a[:, 2], a[:, 1])
    volume = parallelpiped_volume(np.array([a[:, 0], a[:, 2], a[:, 1]]))

    assert np.isclose(triple_product_result, -65.5360), "triple product can be negative"
    assert np.isclose(volume, 65.5360), "volume is always defined as positve"


def test_reciprocal_lattice_vectors():
    r"""
    Test
        \mathbf{a}_i \cdot \mathbf{b}_j = 2 \pi \delta_{ij}

    """

    # fcc with arbitrary lattice constant
    a = 3.2 * np.array([[0., 1., 1.],
                        [1., 0., 1.],
                        [1., 1., 0.]])
    b = reciprocal_lattice_vectors(a)

    assert a.shape == (3, 3)
    assert b.shape == (3, 3)

    a_dot_b = np.transpose(a) @ b
    off_diagonal_indices = np.where(~np.eye(a_dot_b.shape[0], dtype=bool))
    off_diagonal_elements = a_dot_b[off_diagonal_indices]

    assert np.allclose(a_dot_b.diagonal(), 2 * np.pi)
    assert np.allclose(off_diagonal_elements, 0.)
