"""
Functions that operate on lattice vectors 
"""
import numpy as np

from .math_utils import triple_product
import xml.etree.cElementTree as ET


def parse_lattice_vectors(file_path:str) -> np.array:
    """
    Parse the lattice coordinate from the input.xml. Units are in bohr.
    :param file_path:    relative path to the input.xml
    :return lattvec:    matrix that holds the lattice vectors.
    """
    file_name = 'input.xml'
    treeinput = ET.parse(file_path)
    if file_path.split('/')[-1] != file_name:
        file_path = os.path.join(file_path, fileName)

    treeinput = ET.parse(file_path)
    root_input = treeinput.getroot()

    try:
        scale = float(root_input.find("structure").find("crystal").attrib["scale"])
    except Exception:
        scale = 1.0
    lat_vec = [np.array(val.text.split(), dtype=float)*scale for val in root_input.find("structure").find("crystal").findall("basevect")]
    lat_vec = np.array(lat_vec).transpose()
    
    return lat_vec



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


def plane_transformation(rec_lat_vec:np.array, plot_vec:np.array) -> np.array:
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
                           / norm(np.cross(plot_vec[1] - plot_vec[0], plot_vec[2] - plot_vec[0])) ])
    transformation_matrix = np.linalg.inv(plot_vec)
    for v in transformation_matrix:
        v = v / norm(v)
    transformation_matrix = np.transpose(transformation_matrix)

    return transformation_matrix
