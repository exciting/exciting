""" GW eigenvalue parser, processing dict values and returning to an object.
"""
import os
from typing import Optional, List, Union, Dict
import enum
import numpy as np

from excitingtools.exciting_dict_parsers.gw_eigenvalues_parser import parse_evalqp, _file_name, parse_gw_dos
from excitingtools.dataclasses.data_structs import NumberOfStates
from excitingtools.dataclasses.eigenvalues import EigenValues
from excitingtools.dataclasses.density_of_states import DOS


class NitrogenEvalQPColumns(enum.Enum):
    E_KS = 0
    E_HF = 1
    E_GW = 2
    Sx = 3
    Sc = 4
    Vxc = 5
    DE_HF = 6
    DE_GW = 7
    Znk = 8


class OxygenEvalQPColumns(enum.Enum):
    """Columns of `_file_name`, for exciting oxygen
    excluding the state index."""
    E_KS = 0
    E_HF = 1
    E_GW = 2
    Sx = 3
    Re_Sc = 4
    Im_Sc = 5
    Vxc = 6
    DE_HF = 7
    DE_GW = 8
    Znk = 9


# Columns pass as a single enum, or list of enums.
columns_type = Union[enum.Enum, List[enum.Enum]]

# Return an instance of EigenValues, or a dicts of EigenValues if the user asks for
# multiple data columns.
return_type = Union[EigenValues, Dict[enum.Enum, EigenValues]]


def gw_eigenvalue_parser(input_file_path: str, columns: Optional[columns_type] = OxygenEvalQPColumns.E_GW) -> \
        return_type:
    """ High-level Parser for GW eigenvalues file.

    Unpacks the result of dict into a sensible form and returns the data to return_type.

    :param input_file_path: File path (can include or exclude file name).
    :param columns: Optional choice for which data column of energies to return. Default is to return
    GW eigenvalues, assuming exciting oxygen (most recent release).
    :return An instance of EigenValues, or a dict of (key, value) = [EvalQPColumns, EigenValues].
    """
    path, file_name = os.path.split(input_file_path)
    if file_name == _file_name:
        file_path = path
    else:
        file_path = input_file_path
    abs_file_name = os.path.join(file_path, _file_name)

    if not isinstance(columns, list):
        columns = [columns]

    # Parse data
    data: dict = parse_evalqp(abs_file_name)
    state_range_list: List[int] = data.pop('state_range')
    column_enums = data.pop('column_labels')

    n_k = len(data.keys())
    state_range = NumberOfStates(state_range_list[0], state_range_list[1])

    # Check for inconsistency between requested column and columns available (i.e. data produced with
    # different code versions)
    # Because the parsed column data class is dynamically-generated, one has to do it with lengths
    parsed_column_names = {enum.value for enum in column_enums}
    requested_column_names = {enum.value for enum in type(columns[0])}

    if (requested_column_names - parsed_column_names) != set():
        enum_class_name = type(columns[0]).__name__
        raise ValueError(f'The requested data column is indexed according to exciting version {enum_class_name},'
                         f'which is not consistent with the columns of the parsed data.'
                         f' Check that your data was produced with the same code version.')

    # Repackage data
    k_indices = []
    k_points = []
    weights = []
    n_columns = len(type(columns[0]))
    all_eigenvalues = np.empty(shape=(n_k, state_range.n_states, n_columns))

    for ik, k_block in data.items():
        k_indices.append(ik)
        k_points.append(k_block['k_point'])
        weights.append(k_block['weight'])
        all_eigenvalues[ik - 1, :, :] = k_block['energies'][:, :]

    # Return data
    if len(columns) == 1:
        return EigenValues(state_range, k_points, k_indices, all_eigenvalues[:, :, columns[0].value], weights)

    eigen_values = {}
    for column in columns:
        # Note, bad memory access pattern for all_eigenvalues, but would need to transpose in parse_evalqp, then
        # have bad access in the ik loop above
        value = EigenValues(state_range, k_points, k_indices, all_eigenvalues[:, :, column.value], weights)
        eigen_values[column] = value

    return eigen_values


def parse_obj_gw_dos(full_file_name: str) -> DOS:
    """High-level parser for GW DOS files.

    :param full_file_name: Path + file name
    :return: DOS object
    """
    gw_dos_data = parse_gw_dos(full_file_name)
    return DOS(gw_dos_data['energy'], gw_dos_data['dos'])
