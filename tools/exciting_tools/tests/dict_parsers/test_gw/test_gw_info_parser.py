"""Test parsers that comprise the GW_INFO.OUT parser.

One notes that GW_INFO.OUT is plain text, and so any change to the formatting
or text which is used in the parser pattern-matching, will break the parser.
A much better format for this file would be YAML.
"""

import numpy as np

from excitingtools.exciting_dict_parsers.gw_info_parser import parse_gw_info, parse_frequency_grid, \
    parse_correlation_self_energy_params, extract_kpoints, parse_ks_eigenstates, \
    parse_n_q_point_cycles, parse_band_structure_info, parse_mixed_product_params, \
    parse_bare_coulomb_potential_params, parse_gw_timings

# Text files are large, it's easier to store externally.
# Note, these fixtures are used, even if greyed out by the IDE
from . mock_gw_info_out import zro2_gw_info_out_mock, si_2_gw_info_out_mock


def test_parse_correlation_self_energy_params(zro2_gw_info_out_mock):
    """ Test `Correlation self-energy parameters` block of GW_INFO.OUT.
    """
    reference = {'Solution of the QP equation': 0,
                 'Energy alignment': 0,
                 'Analytic continuation method': "PADE - Thiele's reciprocal difference method",
                 'Analytic continuation method citation': "H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179 (1977)",
                 'Scheme to treat singularities': 'Auxiliary function method "mpb"',
                 'Scheme to treat singularities citation': 'S. Massidda, M. Posternak, and A. Baldereschi, PRB 48, 5058 (1993)'
                 }

    output = parse_correlation_self_energy_params(zro2_gw_info_out_mock.string)
    assert reference == output, "parse_correlation_self_energy_params dictionary not consistent with reference"


def test_parse_mixed_product_params(zro2_gw_info_out_mock):
    """ Test `Mixed product basis parameters` block of GW_INFO.OUT.
    """
    ref = {
        'MT Angular momentum cutoff': 4,
        'MT Linear dependence tolerance factor': 0.001,
        'Plane wave cutoff (in units of Gkmax)': 1.0
    }

    output = parse_mixed_product_params(zro2_gw_info_out_mock.string)

    assert output == ref, "Expect parsed mixed product basis parameters to equal the reference"


def test_parse_bare_coulomb_potential_params(zro2_gw_info_out_mock):
    """ Test `Bare Coulomb potential parameters` block of GW_INFO.OUT.
    """
    ref = {
        'Plane wave cutoff (in units of Gkmax*gmb)': 2.0,
        'Error tolerance for structure constants': 1e-16,
        'MB tolerance factor': 0.1
    }

    output = parse_bare_coulomb_potential_params(zro2_gw_info_out_mock.string)

    assert output == ref, "Expect parsed bare Coulomb parameters to match the reference"


def test_parse_frequency_grid(zro2_gw_info_out_mock):
    """ Test parsing `frequency grid` block of GW_INFO.OUT.
    """
    n_points = 32
    f_grid = parse_frequency_grid(zro2_gw_info_out_mock.string, n_points)

    ref_frequencies = np.array([5.2995325042E-03, 2.7712488463E-02, 6.7184398806E-02, 0.1222977958, 0.1910618778,
                                0.2709916112, 0.3591982246, 0.4524937451, 0.5475062549, 0.6408017754,
                                0.7290083888, 0.8089381222, 0.8777022042, 0.9328156012, 0.9722875115,
                                0.9947004675, 1.005327767, 1.028502360, 1.072023237, 1.139338599,
                                1.236188495, 1.371726328, 1.560544990, 1.826463152, 2.209975300,
                                2.783978125, 3.690151129, 5.233906478, 8.176762249, 14.88440795,
                                36.08481430, 188.6958895])

    ref_weights = np.array([1.3576229706E-02, 3.1126761969E-02, 4.7579255841E-02, 6.2314485628E-02, 7.4797994408E-02,
                            8.4578259698E-02, 9.1301707522E-02, 9.4725305228E-02, 9.4725305228E-02, 9.1301707522E-02,
                            8.4578259698E-02, 7.4797994408E-02, 6.2314485628E-02, 4.7579255841E-02, 3.1126761969E-02,
                            1.3576229706E-02, 1.3721277051E-02, 3.2926421206E-02, 5.4679689940E-02, 8.0889962951E-02,
                            0.1143034524, 0.1591452545, 0.2223471090, 0.3160005534, 0.4626375217,
                            0.7076370069, 1.151720377, 2.048999581, 4.166311667, 10.54097479,
                            40.53058703, 483.3971183])

    assert len(ref_frequencies) == 32, "Require 32 reference frequency points"
    assert len(ref_weights) == 32, "Require 32 reference weights"

    assert np.allclose(f_grid[0, :],
                       ref_frequencies), "Frequency points parsed from gw_info_out disagree with reference"
    assert np.allclose(f_grid[1, :], ref_weights), "Weights parsed from gw_info_out disagree with reference"


def test_parse_ks_eigenstates(zro2_gw_info_out_mock):
    """ Test parsing ` Kohn-Sham eigenstates summary` block from GW_INFO.OUT.
    """
    ref = {
        'Maximum number of LAPW states': 847,
        'Minimal number of LAPW states': 838,
        'Number of states used in GW - total KS': 838,
        'Number of states used in GW - occupied': 21,
        'Number of states used in GW - unoccupied': 2000,
        'Number of states used in GW - dielectric function': 838,
        'Number of states used in GW - self energy': 838,
        'Energy of the highest unoccupied state': 1030.791933,
        'Number of valence electrons': 42,
        'Number of valence electrons treated in GW': 42
    }

    output = parse_ks_eigenstates(zro2_gw_info_out_mock.string)

    assert output == ref, "Parsed KS eigenstates settings not consistent with reference"


def test_parse_n_q_point_cycles(zro2_gw_info_out_mock):
    max_q = parse_n_q_point_cycles(zro2_gw_info_out_mock.string)
    assert max_q == 2, "Two q cycles expected from the reference data"


def test_extract_kpoint(zro2_gw_info_out_mock):
    ref = {'VBM': {'k_point': [0.0, 0.5, 0.5], 'ik': 3},
           'CBm': {'k_point': [0.0, 0.0, 0.0], 'ik': 1}
           }

    output = extract_kpoints(zro2_gw_info_out_mock.string)

    assert output == ref, "Expect extracted VBM and CBm k-points to match reference"


def test_parse_band_structure_info(zro2_gw_info_out_mock):
    """ Test parsing the `Kohn-Sham band structure` block of GW_INFO.OUT
    """
    ks_ref = {
        'Fermi energy': 0.0,
        'Energy range': [-14.6863, 1030.7919],
        'Band index of VBM': 21,
        'Band index of CBm': 22,
        'Indirect BandGap (eV)': 3.3206,
        'Direct Bandgap at k(VBM) (eV)': 3.7482,
        'Direct Bandgap at k(CBm) (eV)': 3.8653,
        'VBM': {'k_point': [0.0, 0.5, 0.5], 'ik': 3},
        'CBm': {'k_point': [0.0, 0.0, 0.0], 'ik': 1}
    }
    ks_output = parse_band_structure_info(zro2_gw_info_out_mock.string, 'ks')
    assert ks_output == ks_ref, "Expect parsed KS band structure info to match the reference"


def test_parse_g0w0_band_structure_info(si_2_gw_info_out_mock, zro2_gw_info_out_mock):
    """ Test parsing the `G0W0 band structure` block of GW_INFO.OUT
    """
    # Direct gap
    g0w0_ref = {
        'Fermi energy': 0.0176,
        'Energy range': [-0.4799, 0.5045],
        'Band index of VBM': 4,
        'Band index of CBm': 5,
        'Direct BandGap (eV)': 3.2457,
        'VBM': {'k_point': [0.0, 0.0, 0.0], 'ik': 1},
        'CBm': {'k_point': [0.0, 0.0, 0.0], 'ik': 1}
    }
    gw_output = parse_band_structure_info(si_2_gw_info_out_mock.string, 'gw')
    assert gw_output == g0w0_ref, \
        "Expect parsed G0W0 band structure info to match the reference for direct gap "


    # Indirect gap
    g0w0_ref = {
        'Fermi energy': -0.0054,
        'Energy range': [-16.2632, 1031.409],
        'Band index of VBM': 21,
        'Band index of CBm': 22,
        'Indirect BandGap (eV)': 5.392,
        'Direct Bandgap at k(VBM) (eV)': 5.5472,
        'Direct Bandgap at k(CBm) (eV)': 5.9646,
        'VBM': {'k_point': [0.0, 0.5, 0.5], 'ik': 3},
        'CBm': {'k_point': [0.0, 0.0, 0.0], 'ik': 1}
    }
    gw_output = parse_band_structure_info(zro2_gw_info_out_mock.string, 'gw')
    assert gw_output == g0w0_ref, \
        "Expect parsed G0W0 band structure info to match the reference for indirect gap"


def test_parse_gw_info(zro2_gw_info_out_mock):
    """ Test parsing of the whole GW_INFO.OUT
    """

    # Reference, without frequencies_weights
    ref = {'correlation_self_energy_parameters': {'Solution of the QP equation': 0,
                                                  'Energy alignment': 0,
                                                  'Analytic continuation method': "PADE - Thiele's reciprocal difference method",
                                                  'Analytic continuation method citation':
                                                      'H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179 (1977)',
                                                  'Scheme to treat singularities': 'Auxiliary function method "mpb"',
                                                  'Scheme to treat singularities citation':
                                                      'S. Massidda, M. Posternak, and A. Baldereschi, PRB 48, 5058 (1993)'},
           'mixed_product_basis_parameters': {'MT Angular momentum cutoff': 4,
                                              'MT Linear dependence tolerance factor': 0.001,
                                              'Plane wave cutoff (in units of Gkmax)': 1.0},
           'bare_coulomb_potential_parameters': {'Plane wave cutoff (in units of Gkmax*gmb)': 2.0,
                                                 'Error tolerance for structure constants': 1e-16,
                                                 'MB tolerance factor': 0.1},
           'screened_coulomb_potential': 'Full-frequency Random-Phase Approximation',
           'core_electrons_treatment': 'all - Core states are included in all calculations',
           'qp_interval': [1, 2000],
           'n_empty': 2000,
           'q_grid': [2, 2, 2],
           'mixed_product_wf_info': {'Maximal number of MT wavefunctions per atom': 1069,
                                     'Total number of MT wavefunctions': 2733,
                                     'Maximal number of PW wavefunctions': 468,
                                     'Total number of mixed-product wavefunctions': 3201},
           'frequency_grid': {'Type: < fgrid >': 'gauleg2',
                              'Frequency axis: < fconv >': 'imfreq',
                              'Number of frequencies: < nomeg >': 32,
                              'Cutoff frequency: < freqmax >': 1.0},
           'ks_eigenstates_summary': {'Maximum number of LAPW states': 847,
                                      'Minimal number of LAPW states': 838,
                                      'Number of states used in GW - total KS': 838,
                                      'Number of states used in GW - occupied': 21,
                                      'Number of states used in GW - unoccupied': 2000,
                                      'Number of states used in GW - dielectric function': 838,
                                      'Number of states used in GW - self energy': 838,
                                      'Energy of the highest unoccupied state': 1030.791933,
                                      'Number of valence electrons': 42,
                                      'Number of valence electrons treated in GW': 42},
           'ks_band_structure_summary': {'Fermi energy': 0.0,
                                         'Energy range': [-14.6863, 1030.7919],
                                         'Band index of VBM': 21, 'Band index of CBm': 22,
                                         'Indirect BandGap (eV)': 3.3206,
                                         'Direct Bandgap at k(VBM) (eV)': 3.7482,
                                         'Direct Bandgap at k(CBm) (eV)': 3.8653,
                                         'VBM': {'k_point': [0.0, 0.5, 0.5], 'ik': 3},
                                         'CBm': {'k_point': [0.0, 0.0, 0.0], 'ik': 1}},
           'n_q_cycles': 2,
           'g0w0_band_structure_summary': {'Fermi energy': -0.0054,
                                           'Energy range': [-16.2632, 1031.409],
                                           'Band index of VBM': 21,
                                           'Band index of CBm': 22,
                                           'Indirect BandGap (eV)': 5.392,
                                           'Direct Bandgap at k(VBM) (eV)': 5.5472,
                                           'Direct Bandgap at k(CBm) (eV)': 5.9646,
                                           'VBM': {'k_point': [0.0, 0.5, 0.5], 'ik': 3},
                                           'CBm': {'k_point': [0.0, 0.0, 0.0], 'ik': 1}}
           }

    output = parse_gw_info(zro2_gw_info_out_mock.file)

    # frequencies and weights tested in separate unit test
    f_w = output['frequency_grid'].pop('frequencies_weights')

    assert output == ref, "Output from parse_gw_info does not agree with reference dictionary"


def test_parse_gw_timings(zro2_gw_info_out_mock):
    """ Test parsing `GW timing info` block in GW_INFO.OUT
    """

    ref = {
        'Initialization': {
            'Initialization': 15.46,
            'init_scf': 8.38,
            'init_kpt': 0.04,
            'init_eval': 0.02,
            'init_freq': 0.0,
            'init_mb': 6.76
        },
        'Subroutines': {
            'Subroutines': None,
            'calcpmat': 5.12,
            'calcbarcmb': 5.65,
            'BZ integration weights': 18.35
        },
        'Dielectric function': {
            'Dielectric function': 422.09,
            'head': 0.3,
            'wings': 70.14,
            'body (not timed)': 0.0,
            'inversion': 1.72
        },
        'WF products expansion': {
            'WF products expansion': 2525.59,
            'diagsgi': 0.24,
            'calcmpwipw': 0.04,
            'calcmicm': 2.63,
            'calcminc': 2.1,
            'calcminm': 2520.58
        },
        'Self-energy': {'Self-energy': 7069.46, 'calcselfx': 43.33, 'calcselfc': 7026.14},
        'calcvxcnn': {'calcvxcnn': 27.52},
        'input/output': {'input/output': 0.0},
        'Total': {'Total': 7555.78}
    }

    assert parse_gw_timings(zro2_gw_info_out_mock.string) == ref, "Parsed timings do not agree with reference dictionary"
