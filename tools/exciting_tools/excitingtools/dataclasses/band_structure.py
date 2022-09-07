""" Band structure class.
"""
from typing import Tuple, List, Optional
import numpy as np


class BandData:
    ticks_and_labels = Tuple[np.ndarray, List[str]]
    vertex_keys = ['distance', 'label', 'coord']

    def __init__(self,
                 bands: np.ndarray,
                 k_points: Optional[np.ndarray] = None,
                 flattened_k_points: Optional[np.ndarray] = None,
                 vertices: Optional[List[dict]] = None):
        """ Initialise BandData.

        :param bands: Band energies with shape (n_k_points, n_bands)
        :param: k_points: k-points at which the band energies have been computed.
        :param: flattened_k_points: Flattened k-points along which one can plot a band structure.
        :param: vertices: exciting output containing high-symmetry points and symbols along the
        flattened k-path, as defined in exciting.
        """
        self.bands = bands
        self.n_k_points, self.n_bands = self.bands.shape
        self.k_points = k_points
        self.flattened_k_points = flattened_k_points
        self.vertices = vertices
        self.xticks, self.labels = self.band_path()

    def band_path(self) -> ticks_and_labels:
        """ Get an array of points in the k-path that correspond to high symmetry points,
        and a list of their labels.

        vertices expected to have the form
        [{'distance': float, 'label': str, 'coord': [float, float, float]}, ...]
        parsed from exciting's bandstructure.xml file.

        :return: Tuple of NumPy array containing high symmetry points and list containing their labels.
        If vertices is None, return empty containers.
        """
        if self.vertices is None:
            return np.empty(shape=1), []

        assert list(self.vertices[0].keys()) == self.vertex_keys, \
            f'Expect a vertex to have the keys {self.vertex_keys}'

        vertices = [self.vertices[0]["distance"]]
        labels = [self.vertices[0]["label"]]

        for i in range(1, len(self.vertices)):
            vertex = self.vertices[i]["distance"]
            label = self.vertices[i]["label"]

            # Handle discontinuities in the band path
            if np.isclose(vertex, vertices[-1]):
                vertices.pop()
                label = labels.pop() + ',' + label

            vertices.append(vertex)
            labels.append(label)

        # Replace for plotting purposes
        unicode_gamma = '\u0393'
        for label in ['Gamma', 'gamma', 'G']:
            labels = list(map(lambda x: x.replace(label, unicode_gamma), labels))

        return np.asarray(vertices), labels

    # TODO(Mara) Issue #135. Implement method
    def get_valence_band_maximum(self):
        # Note, difference in indexing between python and fortran
        # self.i_vbm =
        raise NotImplementedError('get_vbm not implemented')

    # TODO(Mara) Issue #135. Implement method
    def get_conduction_band_minimum(self):
        # self.i_cbm =
        raise NotImplementedError('get_cbm not implemented')

    # TODO(Mara) Issue #135. Implement method
    def get_fundamental_band_gap(self):
        # Get k index for for i_vbm
        # Get k index for i_cbm
        # return bands[self.i_cbm - 1, ik_cbm] - bands[self.i_vbm - 1, ik_vbm]
        raise NotImplementedError('get_fundamental_band_gap not implemented')

    # TODO(Mara) Issue #135. Implement method
    def get_band_gap(self, k_valence, k_conduction):
        # Get k index for k_valence
        # Get k index for k_conduction
        # return bands[self.i_cbm - 1, ik_cbm] - bands[self.i_vbm - 1, ik_vbm]
        raise NotImplementedError('get_band_gap not implemented')
