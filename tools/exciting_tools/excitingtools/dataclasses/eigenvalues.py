""" Eigenvalue class.
"""
import warnings
from typing import List, Union, Optional
import numpy as np

from excitingtools.dataclasses.data_structs import PointIndex, BandIndices, NumberOfStates


class EigenValues:
    point_type = Union[np.ndarray, List[float]]
    index_type = Union[np.ndarray, List[int]]
    # If a k-index is not matched
    NO_MATCH = -1

    def __init__(self,
                 state_range: NumberOfStates,
                 k_points: point_type,
                 k_indices: index_type,
                 all_eigenvalues: np.ndarray,
                 weights=None):
        self.state_range = state_range
        self.k_points = k_points
        self.k_indices = k_indices
        self.all_eigenvalues = all_eigenvalues
        self.weights = weights
        if all_eigenvalues.shape != (len(self.k_points), self.state_range.n_states):
            raise ValueError('Shape of all_eigenvalues does not match (n_k, n_states)')

    def get_array_index(self, i_state: int):
        """ Given the state index, get the corresponding index in the eigenvalue array.

        :param i_state: State index using fortran indexing
        :return array index, using zero-indexing
        """
        assert i_state > 0, "state indexing starts at 1"
        return i_state - self.state_range.first_state

    def get_k_point(self, ik: int):
        """ Get the k-point associated with ik index.

        :param ik: k-point index.
        :return k-point in fractional coordinates.
        """
        assert ik > 0, "ik indexing starts at 1"
        return self.k_points[ik - 1]

    def get_index(self, k_point, verbose=False) -> int:
        """ Find the corresponding index of a k-point.

        If no k-point is found, NO_MATCH is returned.

        :param k_point: k-point in fractional coordinates.
        :param verbose: Print warning, if no k-point found.
        :return ik: Corresponding index w.r.t. exciting.
        """
        diff = np.empty(shape=len(self.k_points))
        k_point = np.asarray(k_point)
        for i, point in enumerate(self.k_points):
            diff[i] = np.linalg.norm(np.asarray(point) - k_point)
        indices = np.argwhere(diff < 1.e-8)

        if len(indices) == 0:
            if verbose:
                warnings.warn(f'{k_point} not present in list of k-points')
            return self.NO_MATCH

        indices = indices[0]

        if len(indices) > 1:
            raise ValueError(f'Found degenerate k-points at {indices}')

        ik = indices[0] + 1
        return ik

    def get_k_points(self) -> List[PointIndex]:
        """K-points and their indices

        Note. Fortran indexing for ik always input and output.
        """
        k_points_and_indices = []
        for ik in self.k_indices:
            k_point = self.get_k_point(ik)
            k_points_and_indices.append(PointIndex(k_point, ik))
        return k_points_and_indices

    def get_eigenvalues(self, ik: Optional[int] = None, k_point=None) -> np.ndarray:
        """Return eigenvalues for column="data" w.r.t. k-point.

        :param ik: k-point index.
        :param k_point: k-point, used to find the k-index if ik is not passed.
        :return eigenvalues at ik k-point.
        """
        if ik is None:
            if k_point is None:
                raise ValueError('Must provide either k-index or k-point')
            ik = self.get_index(k_point)

            if ik == self.NO_MATCH:
                return np.empty(shape=0)

        return self.all_eigenvalues[ik - 1, :]

    def band_gap(self, band_indices: BandIndices, k_points=None, k_indices=None) -> Union[float, ValueError]:
        """ Get a band gap for two k-points in the valence band top
        and conduction band bottom, respectively.

        TODO(Alex) Consider adding objects to hold k_points and k_indices.
        However, one can retain support for tuples/lists.

        :param band_indices: Band indices
        :param k_points: k-points for the valence and conduction bands, in that order.
        :param k_indices: k-indices for the valence and conduction bands, in that order.
        :return Band gap. A ValueError is returned if one or both k-points are not present.
        """
        if k_indices is None:
            if k_points is None:
                raise ValueError('Must pass either k_points or k_indices to band_gap method')
            k_indices = (self.get_index(k_points[0]), self.get_index(k_points[1]))

        if self.NO_MATCH in k_indices:
            indices = np.argwhere(k_indices == self.NO_MATCH)[0]
            err_msg = f"".join(f'Requested k-point {k_points[i]} not present. \n' for i in indices)
            return ValueError(err_msg)

        i_vbm = self.get_array_index(band_indices.VBM)
        i_cbm = self.get_array_index(band_indices.CBm)

        ik_vbm = k_indices[0] - 1
        ik_cbm = k_indices[1] - 1

        return self.all_eigenvalues[ik_cbm, i_cbm] - self.all_eigenvalues[ik_vbm, i_vbm]
