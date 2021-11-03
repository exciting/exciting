"""
Module containing parsers for GW output files:
  * GW_INFO.OUT
  * EVALQP.DAT
  * EFERMI_GW.OUT
  * VXCNN.DAT
  * EPS00_GW.OUT
"""
from xml.etree.ElementTree import ParseError
import numpy as np
from typing import List, Callable
import sys
import copy
import re
from pathlib import Path

from excitingtools.utils import get_new_line_indices

from .regex_parser import parse_value_regex, parse_values_regex
from .simple_parser import match_current_return_line_n, match_current_extract_from_line_n
from .grep_parser import grep
from .parser_utils import generic_parser
from ..utils import can_be_float, convert_to_literal


def parse_correlation_self_energy_params(file_string: str) -> dict:
    """
    Parse correlation self-energy parameters.

    Match the key on one line and return the value from the second, according to some
    specified extraction behaviour:

    Solution of the QP equation:
      0 - perturbative solution
    Energy alignment:
      0 - no alignment

    returns {'Solution of the QP equation': 0, 'Energy alignment': 0}

    :param str file_string: Input string
    :return dict data: Matched data
    """
    keys_extractions = {'Solution of the QP equation': lambda x: int(x.split()[0]),
                        'Energy alignment': lambda x: int(x.split()[0]),
                        'Analytic continuation method' : lambda x: x.strip(),
                        'Scheme to treat singularities': lambda x: x.strip()
                        }
    return match_current_extract_from_line_n(file_string, keys_extractions)


def parse_mixed_product_params(file_string: str) -> dict:
    """
    Parse mixed product basis parameters.

    :param str file_string: Input string
    :return dict data: Matched data
    """
    search_keys = ['Angular momentum cutoff:',
                   'Linear dependence tolerance factor:',
                   'Plane wave cutoff \\(in units of Gkmax\\):'
                   ]

    data = parse_values_regex(file_string, search_keys)

    # Rename keys to give better context
    modified_data = {'MT Angular momentum cutoff': data['Angular momentum cutoff'],
                    'MT Linear dependence tolerance factor': data['Linear dependence tolerance factor'],
                     'Plane wave cutoff (in units of Gkmax)': data['Plane wave cutoff (in units of Gkmax)']
                     }

    return modified_data


def parse_bare_coulomb_potential_params(file_string: str) -> dict:
    """
    Extract bare Coulomb parameter data of the form:

    Bare Coulomb potential parameters:
        Plane wave cutoff (in units of Gkmax*input%gw%MixBasis%gmb): 2.00000000000000
        Error tolerance for structure constants:   1.000000000000000E-016
        Tolerance factor to reduce the MB size based on
        the eigenvectors of the bare Coulomb potential:   0.100000000000000

    and return a dictionary of the form:

    {'Plane wave cutoff (in units of Gkmax*gmb)': 2.0,
     'Error tolerance for structure constants': 1e-16,
      'MB tolerance factor': 0.1
    }

    :param str file_string: Input string
    :return dict data: Matched data
    """

    pw_cutoff = list(parse_value_regex(file_string, 'Plane wave cutoff \\(in units of Gkmax\\*input%gw%MixBasis%gmb\\):').values())
    assert len(pw_cutoff) == 1, "Matched plane wave cutoff is ambiguous - more than one match"
    # Defined with less-verbose key
    data = {'Plane wave cutoff (in units of Gkmax*gmb)': float(pw_cutoff[0])}

    data2 = parse_value_regex(file_string, 'Error tolerance for structure constants:')

    data3 = parse_value_regex(file_string, 'the eigenvectors of the bare Coulomb potential:')

    # Defined with less-verbose key (see docs above for full key, over 2 lines)
    modified_data3 = {'MB tolerance factor': data3.pop('the eigenvectors of the bare Coulomb potential')}

    return {**data, **data2, **modified_data3}


def parse_mixed_product_wf_info(file_string: str) -> dict:
    """
    Parse mixed product wave function information.

    :param str file_string: Input string
    :return dict data: Matched data
    """
    wf_info_keys = ['Maximal number of MT wavefunctions per atom:',
                    'Total number of MT wavefunctions:',
                    'Maximal number of PW wavefunctions:',
                    'Total number of mixed-product wavefunctions:'
                    ]

    return parse_values_regex(file_string, wf_info_keys)


def parse_frequency_grid_info(file_string: str) -> dict:
    """
    Parse frequency grid information.

    :param str file_string: Input string
    :return dict data: Matched data
    """
    fgrid_keys = ['Type: < fgrid >',
                  'Frequency axis: < fconv >',
                  'Number of frequencies: < nomeg >',
                  'Cutoff frequency: < freqmax >'
                  ]

    return parse_values_regex(file_string, fgrid_keys)


def parse_frequency_grid(file_string: str, n_points: int) -> np.ndarray:
    """
    Parse the frequency grid and weights used with GW

    :param str file_string: Input string
    :param int n_points: Number of frequency grid points
    :return np.ndarray grid_and_weights: Frequency grid and weights in grid[0, :] and grid[1, :],
    respectively
    """
    index = file_string.find('frequency list: < #    freqs    weight >')
    frequency_lines = file_string[index:].split('\n')

    grid_and_weights = np.empty(shape=(2, n_points))
    for i in range(0, n_points):
        index, frequency, weight = frequency_lines[i + 1].split()
        grid_and_weights[0:2, i] = np.array([float(frequency), float(weight)])

    return grid_and_weights


def parse_ks_eigenstates(file_string: str) -> dict:
    """
    Parse information on the KS eigenstates used.
    Note, final keys will have the trailing whitespace stripped.

    :param str file_string: Input string
    :return dict data: Matched data
    """

    # Trailing whitespace in some instances required for match
    ks_eigenstates_keys = ['Maximum number of LAPW states:',
                           'Minimal number of LAPW states:',
                           '- total KS',
                           '- occupied',
                           '- unoccupied ',
                           '- dielectric function',
                           '- self energy',
                           'Energy of the highest unoccupied state: ',
                           'Number of valence electrons:',
                           'Number of valence electrons treated in GW: '
                           ]

    data = parse_values_regex(file_string, ks_eigenstates_keys)

    # Prepend keys that lack context, such as '- total KS'
    modified_data = {}
    prepend_str = "Number of states used in GW"

    for key, value in data.items():
        if key[0] == '-':
            new_key = prepend_str + " " + key.rstrip().rstrip(':')
        else:
            new_key = key.rstrip().rstrip(':')
        modified_data[new_key] = value

    return modified_data


def parse_n_q_point_cycles(file_string: str) -> int:
    """
    Get the maximum number of q iterations performed.

    :param str file_string: Input string
    :return int n_q_cycles:  maximum nunber of q iterations performed
    """
    matches = re.findall('\\(task_gw\\): q-point cycle, iq =' + '(.+?)\n', file_string)
    n_q_cycles = max([int(string.strip()) for string in matches])
    return n_q_cycles


def extract_kpoint(file_string: str) -> dict:
    """
    Parse the substring of the form:

     at k =    0.000   0.500   0.500 ik =     3

    returning a dictionary of the form:

     {'VBM': {'k_point': [0.0, 0.5, 0.5], 'ik': 3},
      'CBm': {'k_point': [0.0, 0.5, 0.5], 'ik': 3}
      }

    :param str file_string: Input string
    :return dict k_data: Matched data, of the form documented
    """
    def parse_k_match(file_string: str, key: str) -> dict:
        data = {}

        try:
            match = re.search(key + '(.+?)\n', file_string)
            k_point_and_index = match.group(1).split()
            data = {'k_point': [float(k) for k in k_point_and_index[:3]],
                    'ik': int(k_point_and_index[-1])}

        except AttributeError:
            print("extract_kpoint. Did not find the key", match)

        return data

    k_data = parse_k_match(file_string, 'at k      =    ')

    return {'VBM': k_data, 'CBm': k_data}


def extract_kpoints(file_string: str) -> dict:
    """
    Parse the substring of the form:

     at k(VBM) =    0.000   0.500   0.500 ik =     3
        k(CBm) =    0.000   0.000   0.000 ik =     1

    returning a dictionary of the form:

     {'VBM': {'k_point': [0.0, 0.5, 0.5], 'ik': 3},
      'CBm': {'k_point': [0.0, 0.0, 0.0], 'ik': 1}
      }

    :param str file_string: Input string
    :return dict k_data: Matched data, of the form documented
    """
    def parse_k_match(file_string: str, key: str) -> dict:
        data = {}

        match_key_to_parser_key = {'at k\\(VBM\\) = ': 'VBM',
                                      'k\\(CBm\\) = ': 'CBm'}

        try:
            match = re.search(key + '(.+?)\n', file_string)
            parser_key = match_key_to_parser_key[key].replace('\\', "")
            k_point_and_index = match.group(1).split()
            data[parser_key] = {'k_point': [float(k) for k in k_point_and_index[:3]],
                                'ik': int(k_point_and_index[-1])}

        except AttributeError:
            print("extract_kpoints. Did not find the key", match)

        return data

    k_data = parse_k_match(file_string, 'at k\\(VBM\\) = ')
    k_data2 = parse_k_match(file_string, 'k\\(CBm\\) = ')

    return {** k_data, **k_data2}


def parse_band_structure_info(file_string: str, bs_type: str, ) -> dict:
    """
    Parse KS or GW band structure information.

    This routine assumes that the KS band structure info will ALWAYS appear
    before the GW band structure info.

    Two situations can occur.

    * 1. Indirect bandgap:

     Indirect BandGap (eV):                    3.3206
     at k(VBM) =    0.000   0.500   0.500 ik =     3
        k(CBm) =    0.000   0.000   0.000 ik =     1
     Direct Bandgap at k(VBM) (eV):            3.7482
     Direct Bandgap at k(CBm) (eV):            3.8653

    * 2. Direct bandgap:

      Direct BandGap (eV):                      2.3903
      at k      =    0.000   0.000   0.000 ik =     1

    :param str bs_type: Band structure type to parse. Either 'ks' or 'gw'
    :param str file_string: Input string
    :return dict k_data: Matched data
    """
    if bs_type == 'ks':
        # Parse first instance of each key, exploiting that Kohn-Sham band structure
        # comes before G0W0 band structure. This ASSUMES fixed structure to the file
        pass

    elif bs_type == 'gw':
        # Find G0W0 band structure in the file, then start parsing from there
        gw_header = ' G0W0 band structure '
        index = file_string.find('G0W0 band structure ')
        file_string = file_string[index:]

    else:
        sys.exit("bs_type must be 'ks' or 'gw'")

    # Indirect BandGap may not be present
    band_structure_keys = ['Fermi energy:',
                           'Energy range:',
                           'Band index of VBM:',
                           'Band index of CBm:']

    data = parse_values_regex(file_string, band_structure_keys)

    # Only present if there's an indirect gap
    indirect_gap = parse_value_regex(file_string, 'Indirect BandGap \\(eV\\):', silent_key_error=True)

    if indirect_gap:
        direct_gap_keys = ['Direct Bandgap at k\\(VBM\\) \\(eV\\):',
                           'Direct Bandgap at k\\(CBm\\) \\(eV\\):']
        data.update(indirect_gap)
        data.update(parse_values_regex(file_string, direct_gap_keys))
        k_point_data = extract_kpoints(file_string)

    else:
        data.update(parse_value_regex(file_string, 'Direct BandGap \\(eV\\):'))
        k_point_data = extract_kpoint(file_string)

    return {**data, **k_point_data}


def parse_gw_info(file_string: str) -> dict:
    """
    Parse data from GW_INFO.OUT

    Timings are not parsed

    :param str file_string: Parsed file string
    :return: dict data: dictionary of parsed data
    """
    data = {}
    data['correlation_self_energy_parameters'] = parse_correlation_self_energy_params(file_string)
    data['mixed_product_basis_parameters'] = parse_mixed_product_params(file_string)
    data['bare_coulomb_potential_parameters'] = parse_bare_coulomb_potential_params(file_string)
    data['screened_coulomb_potential'] =  match_current_return_line_n(file_string, 'Screened Coulomb potential:').strip()
    data['core_electrons_treatment'] = match_current_return_line_n(file_string, 'Core electrons treatment:').strip()
    data['qp_interval'] = parse_value_regex(file_string, 'Interval of quasiparticle states \\(ibgw, nbgw\\):')\
        ['Interval of quasiparticle states (ibgw, nbgw)']
    data['n_empty'] = parse_value_regex(file_string, 'Number of empty states \\(GW\\):')['Number of empty states (GW)']
    data['q_grid'] = parse_value_regex(file_string, 'k/q-points grid:')['k/q-points grid']
    data['mixed_product_wf_info'] = parse_mixed_product_wf_info(file_string)
    data['frequency_grid'] = parse_frequency_grid_info(file_string)
    n_freq_points = data['frequency_grid']['Number of frequencies: < nomeg >']
    data['frequency_grid']['frequencies_weights'] = parse_frequency_grid(file_string, n_freq_points)
    data['ks_eigenstates_summary'] = parse_ks_eigenstates(file_string)
    data['ks_band_structure_summary'] = parse_band_structure_info(file_string, 'ks')
    data['n_q_cycles'] = parse_n_q_point_cycles(file_string)
    data['g0w0_band_structure_summary'] = parse_band_structure_info(file_string, 'gw')
    return data


def extract_gw_timings_as_list(file_string: str) -> List[str]:
    """
    Extract GW timing string block as a list.

    Utilises the fact that timings are returned at the end of GW_INFO.OUT

    :param str file_string: Parsed file string
    :param list timings: GW timings, with each element storing a
    line of timings as a string.
    """
    file_list = file_string.split('\n')
    for i, line in enumerate(reversed(file_list)):
        if 'GW timing info (seconds)' in line:
            index = i - 2
            return file_list[-index:]

    return None


def parse_gw_timings(file_string: str) -> dict:
    """
    Parse timings returned by the GW method

    Assumptions:
     * Initial timing value is 'Initialization'
     * Final timing value is 'Total'

    Returns a dictionary with each time subheading storing
    a dictionary of timing breakdowns.

    If the subheading has a time associated with it, that is stored
    in the breakdowns with the same key. If it has no associated timing,
    the value is None.

    For example:

    'Subroutines': {'Subroutines': None,
                    'calcpmat': 5.12,
                    'calcbarcmb': 5.65,
                    'BZ integration weights': 18.35
                    }

    :param str file_string: Parsed file string
    :return: dict data: dictionary of parsed timings
    """
    timings = extract_gw_timings_as_list(file_string)

    # Parse and store nested timings
    data = {}
    component_time = {}
    root_key = 'Initialization'
    initial_key = root_key

    for line in timings:

        key = line.strip().split(':')[0].rstrip()
        if len(key) == 0: continue
        is_root_key = key[0] != '-'

        if is_root_key and (key != initial_key):
            data[root_key] = copy.deepcopy(component_time)
            component_time.clear()
            root_key = key

        time_str = line.strip().split(':')[-1]
        component_time[key.strip('-').strip()] = float(time_str) if can_be_float(time_str) else None

    # Store last value, Total
    data[root_key] = {root_key: float(time_str)}

    # Remove superfluous key from table formatting
    arb_str = '_________________________________________________________'
    del data[arb_str]

    return data


def parse_efermi_gw(name: str) -> dict:
    """
    Parser for EFERMI_GW.OUT
    """
    data = {}
    try:
        data["EFERMI_GW"] = np.genfromtxt(name)
    except:
        raise ParseError
    return data


def k_points_from_evalqp(full_file_name: str) -> dict:
    """
    Get the irreducible k/q-points and weights from 'EVALQP.DAT'

    These *can*  differ from those reported in 'KPOINTS', depending
    on the choice of ngridq in GW

    TODO(Alex) Issue 85. Refactor to accept a string and use regex

    :param str full_file_name: Path + file name
    :return dict k_points: k-points and their weights.
    """
    raw_kpoints = grep("k-point", full_file_name).splitlines()

    k_points = {}
    for ik, line in enumerate(raw_kpoints):
        kx, ky, kz, w = line.split()[-4:]
        k_points[ik + 1] = {'k_point': [float(k) for k in [kx, ky, kz]], 'weight': float(w)}

    return k_points


def n_states_from_evalqp(file_string: str) -> int:
    """
    Get the total number of states used per k-point, from EVALQP.DAT.
    Expect n_states to be fixed per k-point, so simply extract the final
    state index associated with the last k-point. State indexing begins at 1.

    :param str file_string: Input string
    :return int n_states: Number of states per k-point (occupied plus empty)
    """
    lines = file_string.splitlines()

    for line in reversed(lines):
        if len(line.strip()) > 0:
            return int(line.split()[0])

    return None


def parse_evalqp_blocks(full_file_name: str, k_points: dict, n_states: int) -> dict:
    """
    Parse energy information from EVALQP.dat

    The function expects k-points of the form:

      k_points[ik] = {'k_point': k_point, 'weight': weight}

    where the k-index (ik) follows fortran indexing convention, and is expected to be
    contiguous. The function returns parsed data, with each element of shape (n_states, 10).

    The routine exploits the repeating structure EVALQP.dat:
       kpoint k1 k2 k3 weight
       header line
       1
       .
       .
       n_states

       kpoint k1 k2 k3 weight
       header line
       1
       .
       .
       n_states

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
                                max_rows=n_states
                                )
        # Ignore first column (state index)
        data[ik] = block_data[:, 1:]
        skip_lines += n_states + (header_size + blank_line)

    return data


def parse_evalqp(full_file_name: str) -> dict:
    """
    Parse GW output file EVALQP.DAT

    Parse  and return data of the form:
      data[ik] = {'k_point': k_point, 'weight': weight, 'energies': energies}

    energies have the shape (n_states, 10), where the 10 elements are defined as:
    ('E_KS', 'E_HF', 'E_GW', 'sigma_x', 'Re_sigma_c', 'Im_sigma_c', 'V_xc', 'delta_HF', 'delta_GW', 'Znk')

    :param str full_file_name: Path + file name
    :return dict data: Parsed k-points and energies from EVALQP.DAT
    """
    if not Path(full_file_name).is_file():
        sys.exit("File does not exist:", full_file_name)

    n_states = generic_parser(full_file_name, n_states_from_evalqp)
    k_points = k_points_from_evalqp(full_file_name)
    energies = parse_evalqp_blocks(full_file_name, k_points, n_states)
    assert len(k_points) == len(energies), "Should be a set of energies for each k-point"

    # Repackage energies with their respective k-points
    data = {}
    for ik in range(1, len(k_points) + 1):
        data[ik] = {'k_point': k_points[ik]['k_point'],
                    'weight': k_points[ik]['weight'],
                    'energies': energies[ik]
                    }
    return data


def vkl_from_vxc(full_file_name: str) -> dict:
    """
    Extract vkl (k-points?) from VXCNN.DAT

    Each k-point header is defined like:
     ik=   1    vkl=  0.0000  0.0000  0.0000

    K indices are not extracted as they are always contiguous,
    with indexing starting at 1.

    TODO(Alex) Issue 85. Refactor to take a string and use regex

    :param str full_file_name: Path + file name
    :return dict vkl: vkl/k-points
    """
    raw_data = grep("ik", full_file_name).splitlines()

    vkl = {}
    for ik, line in enumerate(raw_data):
        vkl[ik + 1] = [float(k) for k in line.split()[-3:]]

    return vkl


def parse_vxnn_vectors(full_file_name: str, vkl: dict, n_states: int) -> dict:
    """
    Parsed VXC diagonal matrix elements from VXCNN.DAT

    The routine exploits the repeating file structure:

     ik=   1    vkl=  0.0000  0.0000  0.0000
        1       -2.908349       -0.000000
        2       -2.864103        0.000000
        3       -2.864103       -0.000000
        .
        n_states

     ik=   2    vkl=  0.0000  0.0000  0.5000
        1       -2.908349        0.000000
        2       -2.864100       -0.000000
        3       -2.864100        0.000000
        .
        n_states

    :param str full_file_name: Path + file name
    :param dict vkl: Dictionary of vkl (k-points?)
    :param int n_states: Total number of occupied plus empty states
     Note, this is constant per q-point.

    :return dict data: Parsed VXC diagonal matrix elements, per k-point.
    """
    # File formatting
    header_size = 1
    blank_line = 1

    data = {}
    skip_lines = header_size

    # Must iterate lowest to highest, else data won't match k-points
    for ik in range(1, len(vkl) + 1):
        vxc_vector = np.loadtxt(full_file_name,
                                skiprows=skip_lines,
                                max_rows=n_states
                                )
        # Ignore first column (state index)
        data[ik] = vxc_vector[:, 1:]
        skip_lines += n_states + (header_size + blank_line)

    return data


def parse_vxcnn(full_file_name: str) -> dict:
    """
    Parser for VXCNN.DAT

    :param str full_file_name: Path + file name
    :return dict data: Parsed k-points (labelled as vkl) and the diagonal elements of
     Vxc, per k-point.
    """
    if not Path(full_file_name).is_file():
        sys.exit("File does not exist:", full_file_name)

    # Can use same parser for n_states
    n_states = generic_parser(full_file_name, n_states_from_evalqp)
    vkl = vkl_from_vxc(full_file_name)
    v_xc = parse_vxnn_vectors(full_file_name, vkl, n_states)
    assert len(vkl) == len(v_xc), "Should be a vector of Vxc_NN for each k-point"

    # Repackage Vxc vectors with their respective k-points
    data = {}
    for ik in range(1, len(vkl) + 1):
        data[ik] = {'vkl': vkl[ik], 'v_xc_nn': v_xc[ik]}

    return data


def parse_eps00_frequencies(file_string: str) -> dict:
    """
    Parse frequencies from EPS00_GW.OUT

    :param str file_string: Input string
    :return dict frequencies: Dictionary of frequencues {index: frequency}
    """
    initial_header_size = 3
    block_size = 6

    # Start file_lines on first block header
    file_lines = file_string.splitlines()
    file = file_lines[initial_header_size:]
    n_freq = int(len(file) / block_size)

    frequencies = {}
    j = 0
    for i in range(0, n_freq):
        j = i * block_size
        frequencies[i + 1] = float(file[j].split()[-1])

    assert len(frequencies) == n_freq

    return frequencies


def parse_eps00_blocks(file_string: str, n_freq: int) -> dict:
    """
    Parser for epsilon blocks in EPS00_GW.OUT

    File of the form:

    (dielectric tensor, random phase approximation)

    frequency index and value:      1    0.01985507
    real part, imaginary part below
     152.69882148    0.00000000    0.00000000         0.00000000    0.00000000    0.00000000
       0.00000000  152.69882148    0.00000000         0.00000000    0.00000000    0.00000000
       0.00000000    0.00000000  152.69882148         0.00000000    0.00000000    0.00000000

    frequency index and value:      2    0.02771249
     real part, imaginary part below
        8.22189228    0.00000000    0.00000000        -0.00000000    0.00000000    0.00000000
        0.00000000    8.22189228    0.00000000         0.00000000   -0.00000000    0.00000000
        0.00000000    0.00000000    8.22189228         0.00000000    0.00000000   -0.00000000

    :param str file_string: Input string
    :param int n_freq: Number of frequency points
    :return dict data: Dict containing real and imaginary parts of the dielectric function,
     eps00, for each frequency point. Frequency indexing (keys) start at 1.
    """
    line = get_new_line_indices(file_string)
    assert file_string[line[0]: line[1]].isspace(), "First line of EPS00_GW.OUT must be a whiteline"

    initial_header = 3
    offset = initial_header - 1
    # Header lines + eps00
    repeat_block_size = 6

    data = {}
    file_string = file_string.splitlines()[initial_header:]

    def extract_eps_i(line: str) -> tuple:
        eps_i_re = np.array(line.split()[0:3], dtype=float)
        eps_i_img = np.array(line.split()[3:6], dtype=float)
        return eps_i_re, eps_i_img

    for i_freq in range(0, n_freq):
        i = offset + i_freq * repeat_block_size

        eps_x_re, eps_x_img = extract_eps_i(file_string[i])
        eps_y_re, eps_y_img = extract_eps_i(file_string[i + 1])
        eps_z_re, eps_z_img = extract_eps_i(file_string[i + 2])

        data[i_freq + 1] = {'re': np.array([eps_x_re, eps_y_re, eps_z_re]),
                            'img': np.array([eps_x_img, eps_y_img, eps_z_img])
                            }

    return data


def parse_eps00_gw(file_string: str) -> dict:
    """
    Parser frequency grid and epsilon00 from EPS00_GW.OUT

    :param str file_string: Input string
    :return dict data: Dict containing the frequency, and real and imaginary parts
     of the dielectric function, eps00, for each frequency point.
    """
    frequencies = parse_eps00_frequencies(file_string)
    eps00 = parse_eps00_blocks(file_string, len(frequencies))
    assert len(frequencies) == len(eps00), \
        "Expect eps00 to have n frequencies consistent with the frequency grid"

    # Repackage frequency points and eps00 together
    data = {}
    for i_freq in range(1, len(frequencies) + 1):
        data[i_freq] = {'frequency': frequencies[i_freq], 'eps00': eps00[i_freq]}

    return data
