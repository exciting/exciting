""" Band structure class.
"""
from typing import Tuple, List
import numpy as np


class BandData:
    ticks_and_labels = Tuple[np.ndarray, List[str]]
    vertex_keys = ['distance', 'label', 'coord']

    def __init__(self,
                 k_points: np.ndarray,
                 bands: np.ndarray,
                 vertices: List[dict]):
        """
        Initialise with k-points along a path, optionally defined in vertices,
        and band energies, or initialising by parsing the XML file.
        """
        self.k_points = k_points
        self.bands = bands
        self.vertices = vertices
        self.n_k_points, self.n_bands = self.bands.shape

    def band_path(self) -> ticks_and_labels:
        """ Get an array of points in the k-path that correspond to high symmetry points,
        and a list of their labels.

        vertices expected to have the form
        [{'distance': float, 'label': str, 'coord': [float, float, float]}, ...]

        :return: Tuple of NumPy array containing high symmetry points and list containing their labels.
        """
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
        labels = list(map(lambda x: x.replace('Gamma', unicode_gamma), labels))

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
