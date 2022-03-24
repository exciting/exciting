"""
Math utilities 
"""
import numpy as np


def triple_product(a, b, c) -> float:
    """
    Vector triple product, defined as 
      \mathbf{a} \cdot (\mathbf{b} \wedge \mathbf{c})

    :param a: Vector a 
    :param b: Vector b
    :param c: Vector c
    :return triple product
    """
    return np.dot(a, np.cross(b, c))


def unit_vector(x: np.ndarray) -> np.ndarray:
    """ 
    Unit vector of a vector 'x' 
      \mathbf{\hat{x}} =  \frac{\mathbf{x}}{|\mathbf{x}|}

    :param x: Vector x
    :result: Unit vector of x
    """
    return x / np.linalg.norm(x)
