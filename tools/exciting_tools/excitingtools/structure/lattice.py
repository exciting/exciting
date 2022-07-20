"""Functions that operate on a crystal lattice.
"""
import numpy as np

from excitingtools.math.math_utils import triple_product


def check_lattice(lattice: list):
    """Check lattice vector consistency.

    Checks that the user supplies three different lattice vectors.
    One could also check the lattice vector angles.

    :param list lattice: Lattice vectors
    """
    if len(lattice) != 3:
        raise ValueError('lattice argument expected to have 3 elements')

    for i in range(0, 3):
        if len(lattice[i]) != 3:
            raise ValueError(
                f'lattice vector {i} expected to have 3 components. Instead has {len(lattice[i])}'
            )

    for i, j in [(0, 1), (1, 2), (0, 2)]:
        if np.allclose(lattice[i], lattice[j], atol=1.e-6):
            raise ValueError(
                f'lattice vectors {i} and {j} are numerically equivalent')


def check_lattice_vector_norms(lattice: list, tol=1.e-6):
    """ Check the norm of each lattice vector.

    :param list lattice: Lattice vectors
    :param tol: Optional tolerance for what is considered numerically zero.
    """
    norms = np.empty(shape=3)

    for i in range(0, 3):
        norms[i] = np.linalg.norm(lattice[i])

    zero_norms = np.where(norms < tol)[0]
    for i in zero_norms:
        raise ValueError(f'lattice vector {i} has a norm of zero')


def parallelepiped_volume(lattice_vectors: np.ndarray) -> float:
    """Volume of a parallelepiped cell.

    :param np.ndarray lattice_vectors: Lattice vectors, stored column-wise
    :return: float: Cell volume
    """
    return np.abs(
        triple_product(lattice_vectors[:, 0], lattice_vectors[:, 1],
                       lattice_vectors[:, 2]))


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


# TODO(Bene) This needs cleaning up. Missing any documenting maths. Not at all clear what's happening
def plane_transformation(rec_lat_vec: np.array,
                         plot_vec: np.array) -> np.array:
    """
    Take reciprocal lattice vectors and ONS of a plane in rec. lat. coordinates where the first two vectors span the plane and the third is normal to them
    and calculate a matrix that transforms points in the plane to the xy plane in cartesian coordinates.
    input:
    :param rec_lat_vec:            reciprocal lattice vectors
    :param plot_vec:               ONS of the plotting plane
    :return transformation_matrix: matrix that transforms k and spin vectors to the plot plane
    """
    norm = np.linalg.norm
    # transform plot vec in cartesian coordinates
    plot_vec = (rec_lat_vec.dot(plot_vec)).transpose()
    # extend plot vec to an orthogonal system
    plot_vec = np.array([(plot_vec[1] - plot_vec[0]) / norm(plot_vec[1] - plot_vec[0]),
                         (plot_vec[2] - plot_vec[0]) / norm(plot_vec[2] - plot_vec[0]),
                         np.cross(plot_vec[1] - plot_vec[0], plot_vec[2] - plot_vec[0]) \
                         / norm(np.cross(plot_vec[1] - plot_vec[0], plot_vec[2] - plot_vec[0]))])
    transformation_matrix = np.linalg.inv(plot_vec)
    for v in transformation_matrix:
        v = v / norm(v)
    transformation_matrix = np.transpose(transformation_matrix)

    return transformation_matrix
