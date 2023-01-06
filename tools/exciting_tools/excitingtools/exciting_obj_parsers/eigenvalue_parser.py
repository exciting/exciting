""" Eigenvalue parser, processing dict values and returning to an object.
"""

import numpy as np
from excitingtools.exciting_dict_parsers.groundstate_parser import parse_eigval
from excitingtools.dataclasses.eigenvalues import EigenValues
from excitingtools.dataclasses.data_structs import NumberOfStates


def parse_eigenvalues(root) -> EigenValues:
    """ High-level parser for eigenvalues. Calls dictionary parser to parse information from the "eigval.xml" file.

    :param root: XML file name, XML string or ElementTree.Element as input.
    :return: EigenValues object.
    """

    eigval_dict = parse_eigval(root)
    data = eigval_dict['kpt']

    k_points = []
    for point in data.values():
        k_points.append([float(x) for x in point['vkl'].split()])

    k_indices = list(map(int, data.keys()))

    n_k = len(k_points)
    n_e = len(data['1']['state'])
    eigenvalues = np.empty((n_k, n_e))
    occupations = np.empty((n_k, n_e))

    for ik in range(n_k):
        for ie in range(n_e):
            eigenvalues[ik, ie] = data[str(ik + 1)]['state'][str(ie + 1)]['eigenvalue']
            occupations[ik, ie] = data[str(ik + 1)]['state'][str(ie + 1)]['occupancy']

    return EigenValues(state_range=NumberOfStates(1, occupations.shape[1]),
                       k_points=k_points,
                       k_indices=k_indices,
                       all_eigenvalues=eigenvalues,
                       occupations=occupations)
