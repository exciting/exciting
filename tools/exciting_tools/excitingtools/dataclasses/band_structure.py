""" Band structure class. """

from typing import Tuple, List, Optional, Union
from excitingtools.eigenstates.eigenstates import get_k_point_index
import numpy as np


class BandData:
    ticks_and_labels = Tuple[np.ndarray, List[str]]
    vertex_keys = ['distance', 'label', 'coord']

    def __init__(self,
                 bands: np.ndarray,
                 k_points: np.ndarray,
                 e_fermi: float,
                 flattened_k_points: Optional[np.ndarray] = None,
                 vertices: Optional[List[dict]] = None):
        """ Initialise BandData.

        :param bands: Band energies with shape (n_k_points, n_bands).
        :param: k_points: k-points at which the band energies have been computed.
        :param: e_fermi: Fermi level.
        :param: flattened_k_points: Flattened k-points along which one can plot a band structure.
        :param: vertices: exciting output containing high-symmetry points and symbols along the
        flattened k-path, as defined in exciting.
        """
        self.bands = bands
        self.n_k_points, self.n_bands = self.bands.shape
        self.k_points = k_points
        self.e_fermi = e_fermi
        self.flattened_k_points = flattened_k_points
        self.vertices = vertices
        self.xticks, self.labels = self.band_path()
        self.i_vbm, self.i_cbm = self.get_band_edges()

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

    def get_k_point_index(self, k_point: Union[List[float], np.ndarray]) -> int:

        if len(k_point) != 3:
            raise TypeError('Expected type for k-point: list or NumPy array of length 3')

        return get_k_point_index(k_point, self.k_points, verbose=False)

    def get_band_edges(self) -> Tuple[int, int]:
        """Get indices of the valence and conduction bands.

        :return Tuple of indices for the valence and conduction bands. A ValueError is returned if the Fermi level is
        located higher than the energy bands.
        """
        valence_max_energies = np.empty(len(self.k_points))

        for k in range(len(self.k_points)):
            valence_max_energies[k] = np.amax(self.bands[k][self.bands[k] <= self.e_fermi])

        vbm_energy = np.amax(valence_max_energies)
        ik_vbm = np.argmax(valence_max_energies)
        n_occupied = len(np.where(self.bands[ik_vbm, :] <= vbm_energy)[0])
        i_vbm = n_occupied - 1

        if i_vbm + 1 >= self.n_bands:
            raise ValueError(f'Fermi level {self.e_fermi} larger than highest band energy {np.amax(self.bands)}')

        return i_vbm, i_vbm + 1

    def get_valence_band_maximum(self) -> float:
        """Get the value of the valence band maximum.
        """
        return np.amax(self.bands[:, self.i_vbm])

    def get_conduction_band_minimum(self) -> float:
        """Get the value of the conduction band minimum.
        """
        return np.amin(self.bands[:, self.i_cbm])

    def get_fundamental_band_gap(self) -> float:
        """Get the value of the fundamental band gap.

        :return fundamental band gap. The band gap is set to zero if bands cross the Fermi level.
        """
        zero_band_gap = 0.0
        if np.any(self.bands[:, self.i_vbm] > self.e_fermi) or np.any(self.bands[:, self.i_cbm] <= self.e_fermi):
            return zero_band_gap

        ik_v = np.argmax(self.bands[:, self.i_vbm])
        ik_c = np.argmin(self.bands[:, self.i_cbm])

        return self.bands[ik_c, self.i_cbm] - self.bands[ik_v, self.i_vbm]

    def get_band_gap(self, k_valence, k_conduction):
        """ Get the value of the band gap calculated between two given k-points.

        :param k_valence: k-point for the valence band.
        :param k_conduction: k-point for the valence band.
        :return band gap.
        """
        ik_v = self.get_k_point_index(k_valence)
        ik_c = self.get_k_point_index(k_conduction)

        return self.bands[ik_c, self.i_cbm] - self.bands[ik_v, self.i_vbm]
