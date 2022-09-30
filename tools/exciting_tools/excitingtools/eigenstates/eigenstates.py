""" Function for finding the corresponding index of a k-point.
"""
import numpy as np
import warnings


def get_k_point_index(k_point: np.ndarray, k_points_to_search: np.ndarray, verbose=False) -> int:
    """ Find the corresponding index of a k-point.

    If no k-point is found, NO_MATCH is returned.

    :param k_point: k-point in fractional coordinates.
    :param k_points_to_search: Array of k-points in which k_point is searched for
    :param verbose: Print warning, if no k-point found.
    :return ik: Corresponding index w.r.t. exciting.
    """
    NO_MATCH = np.NaN

    diff = np.empty(shape=len(k_points_to_search))
    k_point = np.asarray(k_point)
    for i, point in enumerate(k_points_to_search):
        diff[i] = np.linalg.norm(np.asarray(point) - k_point)
    indices = np.argwhere(diff < 1.e-8)

    if len(indices) == 0:
        if verbose:
            warnings.warn(f'{k_point} not present in list of k-points')
        return NO_MATCH

    if len(indices) > 1:
        raise ValueError(f'Found degenerate k-points at {indices}'.replace('\n', ','))

    ik = indices[0]
    return ik
