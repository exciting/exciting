"""
Module containing parsers for GW output files
"""
from xml.etree.ElementTree import ParseError
import numpy as np
from typing import List, Callable
import sys
import copy
import re

from .regex_parser import parse_value_regex, parse_values_regex
from .simple_parser import match_current_return_line_n, match_current_extract_from_line_n
from ..utils import can_be_float


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
        Plane wave cutoff (in units of Gkmax*input%gw%MixBasis%gmb):
        2.00000000000000
        Error tolerance for structure constants:   1.000000000000000E-016
        Tolerance factor to reduce the MB size based on
        the eigenvectors of the bare Coulomb potential:   0.100000000000000

    and return a dictionary of the form:

    {'Plane wave cutoff (in units of Gkmax*gmb)': 2.0,
     'Error tolerance for structure constants': 1e-16,
      'MB tolerance factor': 0.1
    }
    """

    pw_cutoff = match_current_return_line_n(file_string, 'Plane wave cutoff (in units of Gkmax*').strip()
    # Defined with less-verbose key
    data = {'Plane wave cutoff (in units of Gkmax*gmb)': float(pw_cutoff)}

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

    :param str file_string: Input string
    :return dict data: Matched data
    """

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
            new_key = prepend_str + " " + key
        else:
            new_key = key
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
            print("Did not find the key", match)

        return data

    k_data = parse_k_match(file_string, 'at k\\(VBM\\) = ')
    k_data2 = parse_k_match(file_string, 'k\\(CBm\\) = ')

    return {** k_data, **k_data2}


def parse_band_structure_info(file_string: str, bs_type: str, ) -> dict:
    """
    Parse KS or GW band structure information.

    This routine assumes that the KS band structure info will ALWAYS appear
    before the GW band structure info.

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

    band_structure_keys = ['Fermi energy:',
                           'Energy range:',
                           'Band index of VBM:',
                           'Band index of CBm:',
                           'Indirect BandGap \\(eV\\):',
                           'Direct Bandgap at k\\(VBM\\) \\(eV\\):',
                           'Direct Bandgap at k\\(CBm\\) \\(eV\\):'
                           ]

    data = parse_values_regex(file_string, band_structure_keys)
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
    data['screened_coulomb_potential'] =  match_current_return_line_n(file_string, 'Screened Coulomb potential:')
    data['core_electrons_treatment'] = match_current_return_line_n(file_string, 'Core electrons treatment:')
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


def parse_evalqp(name: str) -> dict:
    """
    Parser for EVALQP.DAT
    """
    try:
        data = np.genfromtxt(name, skip_header=2)
    except:
        raise ParseError
    out = {"0": list(data[:, 0]), "1": list(data[:, 1]), "2": list(data[:, 2]), "3": list(data[:, 3]),
           "4": list(data[:, 4]), "5": list(data[:, 5]), "6": list(data[:, 6]), "7": list(data[:, 7]),
           "8": list(data[:, 8]), "9": list(data[:, 9]), "10": list(data[:, 10])}

    return out


def parse_vxcnn(name: str) -> dict:
    """
    Parser for VXCNN.DAT
    """
    try:
        data = np.genfromtxt(name, skip_header=1)
    except:
        raise ParseError
    out = {"0": list(data[:, 0]), "1": list(data[:, 1]), "2": list(data[:, 2])}

    return out


def parse_eps00_gw(name: str) -> dict:
    """
    Parser for EPS00_GW.OUT
    """
    try:
        data = np.loadtxt(name, skiprows=3, comments=["r", "f"])
    except:
        raise ParseError
    out = {"0": list(data[:, 0]), "1": list(data[:, 1]), "2": list(data[:, 2]), "3": list(data[:, 3]),
           "4": list(data[:, 4]), "5": list(data[:, 5])}

    return out
