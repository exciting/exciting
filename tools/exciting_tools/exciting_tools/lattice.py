"""
Functions that operate on lattice vectors 
"""
import numpy as np

from .math_utils import triple_product


def parallelpiped_volume(lattice_vectors: np.ndarray) -> float:
    """
    Volume of a parallelpiped cell, defined as the triple product.

    :param np.ndarray lattice_vectors: Lattice vectors, stored column-wise
    :return: float: Cell volume
    """
    return np.abs(triple_product(lattice_vectors[:, 0], lattice_vectors[:, 1], lattice_vectors[:, 2]))


def reciprocal_lattice_vectors(a: np.ndarray) -> np.ndarray:
    r"""
    Get the reciprocal lattice vectors of real-space lattice vectors \{\mathbf{a}\}:

      \mathbf{b}_0 = 2 \pi \frac{\mathbf{a}_1 \wedge \mathbf{a}_2} {\mathbf{a}_0 \cdot (\mathbf{a}_1 \wedge \mathbf{a}_2)}
      \mathbf{b}_1 = 2 \pi \frac{\mathbf{a}_2 \wedge \mathbf{a}_3} {\mathbf{a}_0 \cdot (\mathbf{a}_1 \wedge \mathbf{a}_2)}
      \mathbf{b}_2 = 2 \pi \frac{\mathbf{a}_0 \wedge \mathbf{a}_1} {\mathbf{a}_0 \cdot (\mathbf{a}_1 \wedge \mathbf{a}_2)}

    :param np.ndarray a: Lattice vectors, stored column-wise
    :return: np.ndarray b: Reciprocal lattice vectors, stored column-wise
    """
    volume = triple_product(a[:, 0], a[:, 1], a[:, 2])
    b = np.empty(shape=(3, 3))
    b[:, 0] = 2 * np.pi * np.cross(a[:, 1], a[:, 2]) / volume
    b[:, 1] = 2 * np.pi * np.cross(a[:, 2], a[:, 0]) / volume
    b[:, 2] = 2 * np.pi * np.cross(a[:, 0], a[:, 1]) / volume
    return b
