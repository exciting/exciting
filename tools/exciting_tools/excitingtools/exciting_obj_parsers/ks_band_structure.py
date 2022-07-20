""" KS (ground state) Band structure Parser, returning to an object.
"""
from typing import Dict
import numpy as np

from excitingtools.dataclasses.band_structure import BandData
from excitingtools.parser_utils.parser_decorators import xml_root


@xml_root
def parse_band_structure_to_arrays(root) -> BandData:
    """ Parse KS band structure from bandstructure.xml.

    :param root: Band structure XML file name, XML string or ElementTree.Element as input.
    :return: BandData object, containing the discrete points sampling the k-path, the band energies,
    and the vertex information (high symmetry points along the band path).
    """
    # Split band structure file contents: title, bands and vertices
    bs_xml: Dict[str, list] = {'title': [], 'band': [], 'vertex': []}
    for item in list(root):
        try:
            bs_xml[item.tag].append(item)
        except KeyError:
            raise KeyError(f'Element tag {item.tag} requires implementing in band structure parser')

    n_bands = len(bs_xml['band'])
    first_band = bs_xml['band'][0]
    n_kpts = len(list(first_band))

    # Same set of flattened k-points, per band - so parse once
    k_points_along_band = np.array([point.get('distance') for point in list(first_band)], dtype=float)

    # Read E(k), per band
    band_energies = np.empty(shape=(n_kpts, n_bands))
    for ib, band in enumerate(bs_xml['band']):
        for ik, point in enumerate(list(band)):
            band_energies[ik, ib] = point.get('eval')

    vertices = []
    for element in bs_xml['vertex']:
        vertices.append({'distance': float(element.get('distance')),
                         'label': element.get('label'),
                         'coord':  [float(x) for x in element.get('coord').split()]})

    return BandData(k_points_along_band, band_energies, vertices)
