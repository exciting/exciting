""" Parse GW's EVALQP.DAT file into dicts.
Also seemed like a logical place to add the Fermi-level parser.
"""
import enum
from itertools import count
import numpy as np
import re
import os

from excitingtools.dataclasses.data_structs import NumberOfStates


# GW eigenvalue file name
_file_name = 'EVALQP.DAT'


def parse_efermi_gw(name: str) -> dict:
    """ Parser for EFERMI_GW.OUT.

    :param name: File name
    :return data: GW Fermi level.
    """
    data = {}
    try:
        data["EFERMI_GW"] = np.genfromtxt(name)
    except Exception:
        raise Exception("Numpy parsing error")
    return data


def k_points_from_evalqp(file_string: str) -> dict:
    """Get the irreducible k/q-points and weights from 'EVALQP.DAT'

    These *can*  differ from those reported in 'KPOINTS', depending
    on the choice of ngridq in GW

    :param file_string: File string.
    :return k_points: k-points and their weights.
    """
    k_points_raw: list = re.findall(r'\s*k-point .*$', file_string, flags=re.MULTILINE)
    k_points = {}
    for ik, line in enumerate(k_points_raw):
        kx, ky, kz, w = line.split()[-4:]
        k_points[ik + 1] = {
            'k_point': [float(k) for k in [kx, ky, kz]],
            'weight': float(w)
        }
    return k_points


def n_states_from_evalqp(file_string: str) -> NumberOfStates:
    return n_states_from_file(file_string, n_header=2)


def n_states_from_file(file_string: str, n_header: int) -> NumberOfStates:
    """Get the total number of states used per k-point, from EVALQP.DAT.

    Expect n_states to be fixed per k-point, however GW corrections can apply
    over any state interval [ibgw: nbgw], therefore extract the first state
    index and the last state index from the first and last k-points,
    respectively - for optimal parsing efficiency.

    Could also extract from input.xml.

    :param str file_string: Input string.
    :param n_header Number of header lines.
    :return int n_states: Number of states per k-point (occupied plus empty).
    """
    lines = file_string.splitlines()

    # ibgw = first_state
    first_state = int(lines[n_header].split()[0])

    # nbgw = last_state
    last_state = None
    for line in reversed(lines):
        if len(line.strip()) > 0:
            last_state = int(line.split()[0])
            break

    assert last_state is not None, "Index of final state to have GW correct applied not found"

    return NumberOfStates(first_state, last_state)


def parse_evalqp_blocks(full_file_name: str, k_points: dict, n_states: int) -> dict:
    """Parse energy information from EVALQP.dat.

    The function expects k-points of the form:

      k_points[ik] = {'k_point': k_point, 'weight': weight}

    where the k-index (ik) follows fortran indexing convention, and is expected to be
    contiguous. The function returns parsed data, with each element of shape (n_states, 10).

    The routine exploits the repeating structure EVALQP.dat:
       kpoint k1 k2 k3 weight
       header line
       first state (ibgw)
       .
       .
       last state (nbgw)

       kpoint k1 k2 k3 weight
       header line
       first state (ibgw)
       .
       .
       last state (nbgw)

    to parse all energies per k-point.

    :param str full_file_name: Path + file name
    :param dict k_points: Dictionary of k-points
    :param int n_states: Total number of occupied plus empty states
     Note, this is constant per q-point.

    :return dict data: Parsed energies from EVALQP.dat
    """
    # File formatting
    header_size = 2
    blank_line = 1

    data = {}
    skip_lines = header_size

    # Must iterate lowest to highest, else data won't match k-points
    for ik in range(1, len(k_points) + 1):
        block_data = np.loadtxt(full_file_name,
                                skiprows=skip_lines,
                                max_rows=n_states)
        # Ignore first column (state index)
        data[ik] = block_data[:, 1:]
        skip_lines += n_states + (header_size + blank_line)

    return data


def parse_evalqp(full_file_name: str) -> dict:
    """Parse GW output file EVALQP.DAT

    Parse  and return data of the form:
      data[ik] = {'k_point': k_point, 'weight': weight, 'energies': energies}

    NOTE(Alex) Would be good to transpose this: row-major access.

    For oxygen release:
      energies have the shape (n_states, 10), where the 10 elements are defined as:
      ('E_KS', 'E_HF', 'E_GW', 'sigma_x', 'Re_sigma_c', 'Im_sigma_c', 'V_xc', 'delta_HF', 'delta_GW', 'Znk')

    For nitrogen release:
      energies have the shape (n_states, 10), where the 10 elements are defined as:
      ('E_KS', 'E_HF', 'E_GW', 'Sx', 'Sc', 'Vxc', 'DE_HF', 'DE_GW', 'Znk')

    :param str full_file_name: Path + file name
    :return dict data: Parsed k-points and energies from EVALQP.DAT
    """
    try:
        with open(full_file_name) as f:
            file_string = f.read()
    except FileNotFoundError:
        raise FileNotFoundError(f"{full_file_name} does not exist")

    eval_indexing: NumberOfStates = n_states_from_evalqp(file_string)
    k_points = k_points_from_evalqp(file_string)

    # Uses np to parse, hence easier to pass file name here
    energies = parse_evalqp_blocks(full_file_name, k_points, eval_indexing.n_states)
    assert len(k_points) == len(energies), "Should be a set of energies for each k-point"

    data = {'state_range': [eval_indexing.first_state, eval_indexing.last_state],
            'column_labels': parse_column_labels(file_string)}

    # Repackage energies with their respective k-points
    for ik in range(1, len(k_points) + 1):
        data[ik] = {
            'k_point': k_points[ik]['k_point'],
            'weight': k_points[ik]['weight'],
            'energies': energies[ik]
        }
    return data


def parse_column_labels(file_string: str) -> enum.Enum:
    """ Parse the column labels of EVALQP.DAT, which vary between code versions.

    :param str file_string: Input string
    return An enum class with the column labels as attributes, with corresponding
    values starting from 0.
    """
    column_str = file_string.splitlines()[1:2][0]
    column_labels = column_str.split()[1:]
    # zip and count ensure enum indexing starts at 0
    return enum.Enum(value='EvalQPColumns', names=zip(column_labels, count()))


def parse_gw_dos(full_file_name: str) -> dict:
    """Parser for GW DOS files.

    :param full_file_name: Path + file name
    :return dict data: Parsed energies and DOS from GW DOS files
    """
    valid_file_names = ['TDOS.OUT', 'TDOS-QP.OUT']
    path, file_name = os.path.split(full_file_name)
    if file_name not in valid_file_names:
        raise ValueError(f'{file_name} not a valid DOS file name.')

    data = np.genfromtxt(full_file_name)
    return {'energy': data[:, 0], 'dos': data[:, 1]}
