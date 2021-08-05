"""
Test parser for GW_INFO.OUT
"""
import pytest
import numpy as np

from excitingtools.parser.gw_parser import parse_gw_info, parse_frequency_grid, \
    parse_correlation_self_energy_params, extract_kpoint, parse_ks_eigenstates, \
    parse_n_q_point_cycles, parse_band_structure_info, parse_mixed_product_params, \
    parse_bare_coulomb_potential_params, parse_gw_timings


gw_info_out = """
================================================================================
=                            GW input parameters                               =
================================================================================
 
 
 GW taskname:
 
   g0w0 - G0W0 run
 
--------------------------------------------------------------------------------
 
 Frequency integration parameters:
 Number of frequencies:           32
 Cutoff frequency:    1.00000000000000     
 Grid type:
   gauleg2 - Double Gauss-Legendre grid: [0, freqmax] + [freqmax, infty]
 Convolution method:
   imfreq : weights calculated for imaginary frequecies
 
--------------------------------------------------------------------------------
 
 Correlation self-energy parameters:
 Solution of the QP equation:
   0 - perturbative solution
 Energy alignment:
   0 - no alignment
 Analytic continuation method:
  PADE - Thiele's reciprocal difference method (by H. J. Vidberg and J. W. Seren
 ce, J. Low Temp. Phys. 29, 179 (1977))
 Scheme to treat singularities:
  Auxiliary function method by S. Massidda, M. Posternak, and A. Baldereschi, PR
 B 48, 5058 (1993)
 
--------------------------------------------------------------------------------
 
 Mixed product basis parameters:
   MT part:
     Angular momentum cutoff:            4
     Linear dependence tolerance factor:   1.000000000000000E-003
   Interstitial:
     Plane wave cutoff (in units of Gkmax):    1.00000000000000     
 
--------------------------------------------------------------------------------
 
 Bare Coulomb potential parameters:
   Plane wave cutoff (in units of Gkmax*input%gw%MixBasis%gmb): 
   2.00000000000000     
   Error tolerance for structure constants:   1.000000000000000E-016
   Tolerance factor to reduce the MB size based on
   the eigenvectors of the bare Coulomb potential:   0.100000000000000     
 
--------------------------------------------------------------------------------
 
 Screened Coulomb potential:
   Full-frequency Random-Phase Approximation
 
--------------------------------------------------------------------------------
 
 Core electrons treatment:
   all - Core states are included in all calculations
 
--------------------------------------------------------------------------------
 
 Interval of quasiparticle states (ibgw, nbgw):       1   2000
 
 Number of empty states (GW):         2000
 
 k/q-points grid:            2           2           2
 
--------------------------------------------------------------------------------
-                           Mixed product WF info                              -
--------------------------------------------------------------------------------
 
  Maximal number of MT wavefunctions per atom:         1069
  Total number of MT wavefunctions:                    2733
  Maximal number of PW wavefunctions:                   468
  Total number of mixed-product wavefunctions:         3201
 
 
--------------------------------------------------------------------------------
-                               frequency grid                                 -
--------------------------------------------------------------------------------
 
 Type: < fgrid > gauleg2                                 
 Frequency axis: < fconv > imfreq                                  
 Number of frequencies: < nomeg >           32
 Cutoff frequency: < freqmax >    1.00000000000000     
 frequency list: < #    freqs    weight > 
   1  5.2995325042E-03  1.3576229706E-02
   2  2.7712488463E-02  3.1126761969E-02
   3  6.7184398806E-02  4.7579255841E-02
   4  0.1222977958      6.2314485628E-02
   5  0.1910618778      7.4797994408E-02
   6  0.2709916112      8.4578259698E-02
   7  0.3591982246      9.1301707522E-02
   8  0.4524937451      9.4725305228E-02
   9  0.5475062549      9.4725305228E-02
  10  0.6408017754      9.1301707522E-02
  11  0.7290083888      8.4578259698E-02
  12  0.8089381222      7.4797994408E-02
  13  0.8777022042      6.2314485628E-02
  14  0.9328156012      4.7579255841E-02
  15  0.9722875115      3.1126761969E-02
  16  0.9947004675      1.3576229706E-02
  17   1.005327767      1.3721277051E-02
  18   1.028502360      3.2926421206E-02
  19   1.072023237      5.4679689940E-02
  20   1.139338599      8.0889962951E-02
  21   1.236188495      0.1143034524    
  22   1.371726328      0.1591452545    
  23   1.560544990      0.2223471090    
  24   1.826463152      0.3160005534    
  25   2.209975300      0.4626375217    
  26   2.783978125      0.7076370069    
  27   3.690151129       1.151720377    
  28   5.233906478       2.048999581    
  29   8.176762249       4.166311667    
  30   14.88440795       10.54097479    
  31   36.08481430       40.53058703    
  32   188.6958895       483.3971183    
 
 WARNING(init_dft_eigenvalues) nstdf > nstfv !
 
 
--------------------------------------------------------------------------------
-                       Kohn-Sham eigenstates summary                          -
--------------------------------------------------------------------------------
 
 Maximum number of LAPW states:                      847
 Minimal number of LAPW states:                      838
 Number of states used in GW:
     - total KS                                      838
     - occupied                                       21
     - unoccupied                                   2000
     - dielectric function                           838
     - self energy                                   838
 Energy of the highest unoccupied state:     1030.791933
 Number of valence electrons:                         42
 Number of valence electrons treated in GW:           42
 
 WARNING(init_dft_eigenvalues) One uses the maximum number of available states!
 
 
--------------------------------------------------------------------------------
-                          Kohn-Sham band structure                            -
--------------------------------------------------------------------------------
 
 Fermi energy:     0.0000
 Energy range:   -14.6863 1030.7919
 Band index of VBM:  21
 Band index of CBm:  22
 
 Indirect BandGap (eV):                    3.3206
 at k(VBM) =    0.000   0.500   0.500 ik =     3
    k(CBm) =    0.000   0.000   0.000 ik =     1
 Direct Bandgap at k(VBM) (eV):            3.7482
 Direct Bandgap at k(CBm) (eV):            3.8653
 
================================================================================
=                                  GW cycle                                    =
================================================================================
 
 (task_gw): q-point cycle, iq =            1
 (task_gw): q-point cycle, iq =            2
 
--------------------------------------------------------------------------------
-                            G0W0 band structure                               -
--------------------------------------------------------------------------------
 
 Fermi energy:    -0.0054
 Energy range:   -16.2632 1031.4090
 Band index of VBM:  21
 Band index of CBm:  22
 
 Indirect BandGap (eV):                    5.3920
 at k(VBM) =    0.000   0.500   0.500 ik =     3
    k(CBm) =    0.000   0.000   0.000 ik =     1
 Direct Bandgap at k(VBM) (eV):            5.5472
 Direct Bandgap at k(CBm) (eV):            5.9646
 
================================================================================
=                          GW timing info (seconds)                            =
================================================================================
 
 Initialization                             :        15.46
     - init_scf                             :         8.38
     - init_kpt                             :         0.04
     - init_eval                            :         0.02
     - init_freq                            :         0.00
     - init_mb                              :         6.76
 Subroutines                                : 
     - calcpmat                             :         5.12
     - calcbarcmb                           :         5.65
     - BZ integration weights               :        18.35
     Dielectric function                    :       422.09
     - head                                 :         0.30
     - wings                                :        70.14
     - body (not timed)                     :         0.00
     - inversion                            :         1.72
     WF products expansion                  :      2525.59
     - diagsgi                              :         0.24
     - calcmpwipw                           :         0.04
     - calcmicm                             :         2.63
     - calcminc                             :         2.10
     - calcminm                             :      2520.58
     Self-energy                            :      7069.46
     - calcselfx                            :        43.33
     - calcselfc                            :      7026.14
     calcvxcnn                              :        27.52
     input/output                           :         0.00
 _________________________________________________________
 Total                                      :      7555.78
"""


def test_parse_correlation_self_energy_params():

    reference = {'Solution of the QP equation': 0,
                 'Energy alignment': 0,
                 'Analytic continuation method': "PADE - Thiele's reciprocal difference method (by H. J. Vidberg and J. W. Seren",
                 'Scheme to treat singularities': 'Auxiliary function method by S. Massidda, M. Posternak, and A. Baldereschi, PR'}

    output = parse_correlation_self_energy_params(gw_info_out)

    assert reference == output, "parse_correlation_self_energy_params dictionary not consistent with reference"


def test_parse_mixed_product_params():

    ref = {'MT Angular momentum cutoff': 4,
           'MT Linear dependence tolerance factor': 0.001,
           'Plane wave cutoff (in units of Gkmax)': 1.0
           }

    output = parse_mixed_product_params(gw_info_out)

    assert output == ref, "Expect parsed mixed product basis parameters to equal the reference"


def test_parse_bare_coulomb_potential_params():

    ref = {'Plane wave cutoff (in units of Gkmax*gmb)': 2.0,
           'Error tolerance for structure constants': 1e-16,
           'MB tolerance factor': 0.1
           }

    output = parse_bare_coulomb_potential_params(gw_info_out)

    assert output == ref, "Expect parsed bare Coulomb parameters to match the reference"


def test_parse_frequency_grid():
    n_points = 32
    f_grid = parse_frequency_grid(gw_info_out, n_points)

    ref_frequencies = np.array([5.2995325042E-03, 2.7712488463E-02, 6.7184398806E-02, 0.1222977958, 0.1910618778,
                                0.2709916112,     0.3591982246,     0.4524937451,     0.5475062549, 0.6408017754,
                                0.7290083888,     0.8089381222,     0.8777022042,     0.9328156012, 0.9722875115,
                                0.9947004675,      1.005327767,      1.028502360,      1.072023237,  1.139338599,
                                1.236188495,      1.371726328,      1.560544990,      1.826463152,   2.209975300,
                                2.783978125,      3.690151129,      5.233906478,      8.176762249,   14.88440795,
                                36.08481430,      188.6958895])

    ref_weights = np.array([1.3576229706E-02, 3.1126761969E-02, 4.7579255841E-02, 6.2314485628E-02, 7.4797994408E-02,
                            8.4578259698E-02, 9.1301707522E-02, 9.4725305228E-02, 9.4725305228E-02, 9.1301707522E-02,
                            8.4578259698E-02, 7.4797994408E-02, 6.2314485628E-02, 4.7579255841E-02, 3.1126761969E-02,
                            1.3576229706E-02, 1.3721277051E-02, 3.2926421206E-02, 5.4679689940E-02, 8.0889962951E-02,
                            0.1143034524,     0.1591452545,     0.2223471090,     0.3160005534,     0.4626375217,
                            0.7076370069,     1.151720377,      2.048999581,      4.166311667,      10.54097479,
                            40.53058703,      483.3971183])

    assert len(ref_frequencies) == 32, "Require 32 ref frequency points"
    assert len(ref_weights) == 32, "Require 32 ref weights"

    assert np.allclose(f_grid[0, :], ref_frequencies), "Frequency points parsed from gw_info_out disagree with reference"
    assert np.allclose(f_grid[1, :], ref_weights), "Weights parsed from gw_info_out disagree with reference"


def test_parse_ks_eigenstates():

    ref = {'Maximum number of LAPW states': 847,
           'Minimal number of LAPW states': 838,
           'Number of states used in GW - total KS': 838,
           'Number of states used in GW - occupied': 21,
           'Number of states used in GW - unoccupied ': 2000,
           'Number of states used in GW - dielectric function': 838,
           'Number of states used in GW - self energy': 838,
           'Energy of the highest unoccupied state: ': 1030.791933,
           'Number of valence electrons': 42,
           'Number of valence electrons treated in GW: ': 42
           }

    output = parse_ks_eigenstates(gw_info_out)

    assert output == ref, "Parsed KS eigenstates settings not consistent with reference"



def test_parse_n_q_point_cycles():
    max_q = parse_n_q_point_cycles(gw_info_out)
    assert max_q == 2, "Two q cycles expected from the reference data"


def test_extract_kpoint():

    ref = {'VBM': {'k_point': [0.0, 0.5, 0.5], 'ik': 3},
           'CBm': {'k_point': [0.0, 0.0, 0.0], 'ik': 1}
          }

    output = extract_kpoint(gw_info_out)

    assert output == ref, "Expect extracted VBM and CBm k-points to match reference"

# TODO(ALEX UP TO HERE
def test_parse_band_structure_info():

    # See Kohn-Sham band structure in gw_info_out
    ks_ref = {'Fermi energy': 0.0,
              'Energy range': [-14.6863, 1030.7919],
              'Band index of VBM': 21,
              'Band index of CBm': 22,
              'Indirect BandGap (eV)': 3.3206,
              'Direct Bandgap at k(VBM) (eV)': 3.7482,
              'Direct Bandgap at k(CBm) (eV)': 3.8653,
              'VBM': {'k_point': [0.0, 0.5, 0.5], 'ik': 3},
              'CBm': {'k_point': [0.0, 0.0, 0.0], 'ik': 1}
              }

    # See G0W0 band structure in gw_info_out
    g0w0_ref =  {'Fermi energy': -0.0054,
                 'Energy range': [-16.2632, 1031.409],
                 'Band index of VBM': 21,
                 'Band index of CBm': 22,
                 'Indirect BandGap (eV)': 5.392,
                 'Direct Bandgap at k(VBM) (eV)': 5.5472,
                 'Direct Bandgap at k(CBm) (eV)': 5.9646,
                 'VBM': {'k_point': [0.0, 0.5, 0.5], 'ik': 3},
                 'CBm': {'k_point': [0.0, 0.0, 0.0], 'ik': 1}
                 }

    ks_output = parse_band_structure_info(gw_info_out, 'ks')
    gw_output = parse_band_structure_info(gw_info_out, 'gw')

    assert ks_output == ks_ref, "Expect parsed KS band structure info to match the reference"
    assert gw_output == g0w0_ref, "Expect parsed G0W0 band structure info to match the reference"


def test_parse_gw_info():

    # All the parsers combined

    # Reference, without frequencies_weights
    ref = {'correlation_self_energy_parameters': {'Solution of the QP equation': 0,
                                                  'Energy alignment': 0,
                                                  'Analytic continuation method': "PADE - Thiele's reciprocal difference method (by H. J. Vidberg and J. W. Seren",
                                                  'Scheme to treat singularities': 'Auxiliary function method by S. Massidda, M. Posternak, and A. Baldereschi, PR'},
           'mixed_product_basis_parameters': {'MT Angular momentum cutoff': 4,
                                              'MT Linear dependence tolerance factor': 0.001,
                                              'Plane wave cutoff (in units of Gkmax)': 1.0},
           'bare_coulomb_potential_parameters': {'Plane wave cutoff (in units of Gkmax*gmb)': 2.0,
                                                 'Error tolerance for structure constants': 1e-16,
                                                 'MB tolerance factor': 0.1},
           'screened_coulomb_potential': '   Full-frequency Random-Phase Approximation',
           'core_electrons_treatment': '   all - Core states are included in all calculations',
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
                                      'Number of states used in GW - unoccupied ': 2000,
                                      'Number of states used in GW - dielectric function': 838,
                                      'Number of states used in GW - self energy': 838,
                                      'Energy of the highest unoccupied state: ': 1030.791933,
                                      'Number of valence electrons': 42,
                                      'Number of valence electrons treated in GW: ': 42},
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


    output = parse_gw_info(gw_info_out)
    # frequencies and weights tested in separate unit test
    f_w = output['frequency_grid'].pop('frequencies_weights')

    assert output == ref, "Output from parse_gw_info does not agree with reference dictionary"


def test_parse_gw_timings():

    ref = {
        'Initialization': {'Initialization': 15.46,
                           'init_scf': 8.38,
                           'init_kpt': 0.04,
                           'init_eval': 0.02,
                           'init_freq': 0.0,
                           'init_mb': 6.76},

        'Subroutines': {'Subroutines': None,
                        'calcpmat': 5.12,
                        'calcbarcmb': 5.65,
                        'BZ integration weights': 18.35},

        'Dielectric function': {'Dielectric function': 422.09,
                                'head': 0.3,
                                'wings': 70.14,
                                'body (not timed)': 0.0,
                                'inversion': 1.72},

        'WF products expansion': {'WF products expansion': 2525.59,
                                  'diagsgi': 0.24,
                                  'calcmpwipw': 0.04,
                                  'calcmicm': 2.63,
                                  'calcminc': 2.1,
                                  'calcminm': 2520.58},

        'Self-energy': {'Self-energy': 7069.46,
                        'calcselfx': 43.33,
                        'calcselfc': 7026.14},

        'calcvxcnn': {'calcvxcnn': 27.52},
        'input/output': {'input/output': 0.0},
        'Total': {'Total': 7555.78}
           }

    assert parse_gw_timings(gw_info_out) == ref, "Parsed timings do not agree with reference dictionary"
