""" KS (ground state) Band structure Parser, returning to an object.
"""

from excitingtools.dataclasses.band_structure import BandData
from excitingtools.exciting_dict_parsers.properties_parser import parse_band_structure_xml


def parse_band_structure(root) -> BandData:
    """ High-level parser for KS band structure.

    :param root: Band structure XML file name, XML string or ElementTree.Element as input.
    :return: BandData object.
    """
    # TODO(Mara) Add parsing for bandstructure.dat, to get 3D k-points. Then add below
    data = parse_band_structure_xml(root)
    return BandData(data['band_energies'],
                    flattened_k_points=data['k_points_along_band'],
                    vertices=data['vertices'])
