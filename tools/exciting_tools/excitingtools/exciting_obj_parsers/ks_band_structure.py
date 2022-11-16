""" KS (ground state) Band structure Parser, returning to an object.
"""

import os
import numpy as np
from excitingtools.dataclasses.band_structure import BandData
from excitingtools.exciting_dict_parsers.properties_parser import parse_band_structure_xml, parse_band_structure_dat
from pathlib import Path


def parse_band_structure(file_name: str) -> BandData:
    """ High-level parser for KS band structure. Calls dictionary parsers to parse information from both
    "bandstructure.xml" and "bandstructure.dat" files.

    :param str file_name: File name for band structure such as "bandstructure.dat", "bandstructure.xml" or
                          "bandstructure"
    :return: BandData object.
    """
    name = os.path.splitext(file_name)[0]
    data_from_xml = parse_band_structure_xml(name + '.xml')
    data_from_dat = parse_band_structure_dat(name + '.dat')

    e_fermi = float(np.loadtxt(Path(file_name).parent / "EFERMI.OUT"))

    return BandData(bands=data_from_xml['band_energies'],
                    k_points=data_from_dat['k_points'],
                    e_fermi=e_fermi,
                    flattened_k_points=data_from_xml['k_points_along_band'],
                    vertices=data_from_xml['vertices'])
