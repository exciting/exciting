"""
Test parsers for GW output files
"""
import pytest
import numpy as np
import os
import re

from excitingtools.utils import get_new_line_indices

from excitingtools.parser.gw_parser import parse_gw_info, parse_frequency_grid, \
    parse_correlation_self_energy_params, extract_kpoints, parse_ks_eigenstates, \
    parse_n_q_point_cycles, parse_band_structure_info, parse_mixed_product_params, \
    parse_bare_coulomb_potential_params, parse_gw_timings, k_points_from_evalqp, \
    n_states_from_evalqp, parse_evalqp, parse_vxcnn, vkl_from_vxc, parse_eps00_frequencies, \
    parse_eps00_gw

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
  PADE - Thiele's reciprocal difference method (by H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179 (1977))
 Scheme to treat singularities:
  Auxiliary function method by S. Massidda, M. Posternak, and A. Baldereschi, PRB 48, 5058 (1993)
 
--------------------------------------------------------------------------------
 
 Mixed product basis parameters:
   MT part:
     Angular momentum cutoff:            4
     Linear dependence tolerance factor:   1.000000000000000E-003
   Interstitial:
     Plane wave cutoff (in units of Gkmax):    1.00000000000000     
 
--------------------------------------------------------------------------------
 
 Bare Coulomb potential parameters:
   Plane wave cutoff (in units of Gkmax*input%gw%MixBasis%gmb): 2.00000000000000     
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

    reference = {
        'Solution of the QP equation': 0,
        'Energy alignment': 0,
        'Analytic continuation method': "PADE - Thiele's reciprocal difference method (by H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179 (1977))",
        'Scheme to treat singularities': 'Auxiliary function method by S. Massidda, M. Posternak, and A. Baldereschi, PRB 48, 5058 (1993)'
        }

    output = parse_correlation_self_energy_params(gw_info_out)

    assert reference == output, "parse_correlation_self_energy_params dictionary not consistent with reference"


def test_parse_mixed_product_params():

    ref = {
        'MT Angular momentum cutoff': 4,
        'MT Linear dependence tolerance factor': 0.001,
        'Plane wave cutoff (in units of Gkmax)': 1.0
        }

    output = parse_mixed_product_params(gw_info_out)

    assert output == ref, "Expect parsed mixed product basis parameters to equal the reference"


def test_parse_bare_coulomb_potential_params():

    ref = {
        'Plane wave cutoff (in units of Gkmax*gmb)': 2.0,
        'Error tolerance for structure constants': 1e-16,
        'MB tolerance factor': 0.1
        }

    output = parse_bare_coulomb_potential_params(gw_info_out)

    assert output == ref, "Expect parsed bare Coulomb parameters to match the reference"


def test_parse_frequency_grid():
    n_points = 32
    f_grid = parse_frequency_grid(gw_info_out, n_points)

    ref_frequencies = np.array([
        5.2995325042E-03,
        2.7712488463E-02,
        6.7184398806E-02,
        0.1222977958,
        0.1910618778,
        0.2709916112,
        0.3591982246,
        0.4524937451,
        0.5475062549,
        0.6408017754,
        0.7290083888,
        0.8089381222,
        0.8777022042,
        0.9328156012,
        0.9722875115,
        0.9947004675,
        1.005327767,
        1.028502360,
        1.072023237,
        1.139338599,
        1.236188495,
        1.371726328,
        1.560544990,
        1.826463152,
        2.209975300,
        2.783978125,
        3.690151129,
        5.233906478,
        8.176762249,
        14.88440795,
        36.08481430,
        188.6958895
        ])

    ref_weights = np.array([
        1.3576229706E-02,
        3.1126761969E-02,
        4.7579255841E-02,
        6.2314485628E-02,
        7.4797994408E-02,
        8.4578259698E-02,
        9.1301707522E-02,
        9.4725305228E-02,
        9.4725305228E-02,
        9.1301707522E-02,
        8.4578259698E-02,
        7.4797994408E-02,
        6.2314485628E-02,
        4.7579255841E-02,
        3.1126761969E-02,
        1.3576229706E-02,
        1.3721277051E-02,
        3.2926421206E-02,
        5.4679689940E-02,
        8.0889962951E-02,
        0.1143034524,
        0.1591452545,
        0.2223471090,
        0.3160005534,
        0.4626375217,
        0.7076370069,
        1.151720377,
        2.048999581,
        4.166311667,
        10.54097479,
        40.53058703,
        483.3971183
        ])

    assert len(ref_frequencies) == 32, "Require 32 ref frequency points"
    assert len(ref_weights) == 32, "Require 32 ref weights"

    assert np.allclose(f_grid[0, :], ref_frequencies), "Frequency points parsed from gw_info_out disagree with reference"
    assert np.allclose(f_grid[1, :], ref_weights), "Weights parsed from gw_info_out disagree with reference"


def test_parse_ks_eigenstates():

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

    output = parse_ks_eigenstates(gw_info_out)

    assert output == ref, "Parsed KS eigenstates settings not consistent with reference"


def test_parse_n_q_point_cycles():
    max_q = parse_n_q_point_cycles(gw_info_out)
    assert max_q == 2, "Two q cycles expected from the reference data"


def test_extract_kpoint():

    ref = {
        'VBM': {'k_point': [0.0, 0.5, 0.5], 'ik': 3}, 'CBm': {'k_point': [0.0, 0.0, 0.0], 'ik': 1}
        }

    output = extract_kpoints(gw_info_out)

    assert output == ref, "Expect extracted VBM and CBm k-points to match reference"


def test_parse_band_structure_info():

    # See Kohn-Sham band structure in gw_info_out
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

    # See G0W0 band structure in gw_info_out
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

    ks_output = parse_band_structure_info(gw_info_out, 'ks')
    gw_output = parse_band_structure_info(gw_info_out, 'gw')

    assert ks_output == ks_ref, "Expect parsed KS band structure info to match the reference"
    assert gw_output == g0w0_ref, "Expect parsed G0W0 band structure info to match the reference"


def test_parse_gw_info():

    # All the parsers combined

    # Reference, without frequencies_weights
    ref = {
        'correlation_self_energy_parameters': {
            'Solution of the QP equation': 0,
            'Energy alignment': 0,
            'Analytic continuation method': "PADE - Thiele's reciprocal difference method (by H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179 (1977))",
            'Scheme to treat singularities': 'Auxiliary function method by S. Massidda, M. Posternak, and A. Baldereschi, PRB 48, 5058 (1993)'
            },
        'mixed_product_basis_parameters': {
            'MT Angular momentum cutoff': 4,
            'MT Linear dependence tolerance factor': 0.001,
            'Plane wave cutoff (in units of Gkmax)': 1.0
            },
        'bare_coulomb_potential_parameters': {
            'Plane wave cutoff (in units of Gkmax*gmb)': 2.0,
            'Error tolerance for structure constants': 1e-16,
            'MB tolerance factor': 0.1
            },
        'screened_coulomb_potential': 'Full-frequency Random-Phase Approximation',
        'core_electrons_treatment': 'all - Core states are included in all calculations',
        'qp_interval': [1, 2000],
        'n_empty': 2000,
        'q_grid': [2, 2, 2],
        'mixed_product_wf_info': {
            'Maximal number of MT wavefunctions per atom': 1069,
            'Total number of MT wavefunctions': 2733,
            'Maximal number of PW wavefunctions': 468,
            'Total number of mixed-product wavefunctions': 3201
            },
        'frequency_grid': {
            'Type: < fgrid >': 'gauleg2',
            'Frequency axis: < fconv >': 'imfreq',
            'Number of frequencies: < nomeg >': 32,
            'Cutoff frequency: < freqmax >': 1.0
            },
        'ks_eigenstates_summary': {
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
            },
        'ks_band_structure_summary': {
            'Fermi energy': 0.0,
            'Energy range': [-14.6863, 1030.7919],
            'Band index of VBM': 21,
            'Band index of CBm': 22,
            'Indirect BandGap (eV)': 3.3206,
            'Direct Bandgap at k(VBM) (eV)': 3.7482,
            'Direct Bandgap at k(CBm) (eV)': 3.8653,
            'VBM': {'k_point': [0.0, 0.5, 0.5], 'ik': 3},
            'CBm': {'k_point': [0.0, 0.0, 0.0], 'ik': 1}
            },
        'n_q_cycles': 2,
        'g0w0_band_structure_summary': {
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
        }

    output = parse_gw_info(gw_info_out)
    # frequencies and weights tested in separate unit test
    f_w = output['frequency_grid'].pop('frequencies_weights')

    assert output == ref, "Output from parse_gw_info does not agree with reference dictionary"


def test_parse_gw_timings():

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

    assert parse_gw_timings(gw_info_out) == ref, "Parsed timings do not agree with reference dictionary"


def write_file(file_name: str, string: str):
    fid = open(file_name, mode='w')
    fid.write(string)
    fid.close()


# Sample EVALQP.OUT data
evalqp_str = """k-point #     1:    0.000000    0.000000    0.000000    0.125000
 state   E_KS       E_HF       E_GW       Sx         Re(Sc)     Im(Sc)     Vxc        DE_HF      DE_GW      Znk
   1   -14.68627  -16.50308  -16.31014   -4.72516    0.17351    0.00002   -2.90835   -1.81681   -1.62387    0.98817
   2   -11.48914  -12.99338  -12.79204   -4.36834    0.18150    0.00003   -2.86410   -1.50424   -1.30290    0.98500
   3   -11.48914  -12.99338  -12.79204   -4.36834    0.18150    0.00003   -2.86410   -1.50424   -1.30290    0.98500
   4   -11.48914  -12.99338  -12.79203   -4.36834    0.18151    0.00003   -2.86410   -1.50424   -1.30289    0.98500
   5    -6.24399   -7.20716   -7.01655   -3.72742    0.17144    0.00000   -2.76425   -0.96317   -0.77257    0.97579
   6    -6.24399   -7.20716   -7.01652   -3.72742    0.17144   -0.00001   -2.76425   -0.96317   -0.77253    0.97575
   7    -6.24340   -7.20928   -7.01752   -3.73069    0.17268    0.00003   -2.76481   -0.96588   -0.77412    0.97594
   8    -6.24340   -7.20928   -7.01752   -3.73069    0.17268    0.00003   -2.76481   -0.96588   -0.77412    0.97594
   9    -6.24340   -7.20928   -7.01752   -3.73069    0.17268    0.00003   -2.76481   -0.96588   -0.77412    0.97594
  10    -1.82156   -2.31023   -2.07248   -1.52298    0.20180    0.00623   -1.03431   -0.48867   -0.25092    0.87468
  11    -1.00771   -1.29282   -1.07701   -1.21488    0.19790    0.02221   -0.92977   -0.28511   -0.06930    0.79459

k-point #     2:    0.000000    0.000000    0.500000    0.500000
 state   E_KS       E_HF       E_GW       Sx         Re(Sc)     Im(Sc)     Vxc        DE_HF      DE_GW      Znk
   1   -14.68627  -16.50308  -16.31014   -4.72516    0.17351    0.00002   -2.90835   -1.81681   -1.62387    0.98817
   2   -11.48915  -12.99338  -12.79204   -4.36834    0.18150    0.00003   -2.86410   -1.50424   -1.30289    0.98500
   3   -11.48915  -12.99338  -12.79204   -4.36834    0.18151    0.00003   -2.86410   -1.50424   -1.30289    0.98500
   4   -11.48914  -12.99338  -12.79204   -4.36834    0.18150    0.00003   -2.86410   -1.50424   -1.30290    0.98500
   5    -6.24401   -7.20719   -7.01653   -3.72740    0.17145   -0.00001   -2.76423   -0.96318   -0.77252    0.97574
   6    -6.24401   -7.20719   -7.01650   -3.72740    0.17145   -0.00001   -2.76423   -0.96318   -0.77249    0.97570
   7    -6.24345   -7.20933   -7.01756   -3.73065    0.17268    0.00003   -2.76477   -0.96588   -0.77411    0.97593
   8    -6.24344   -7.20932   -7.01754   -3.73066    0.17268    0.00003   -2.76478   -0.96588   -0.77410    0.97593
   9    -6.24344   -7.20932   -7.01755   -3.73066    0.17268    0.00003   -2.76478   -0.96588   -0.77411    0.97593
  10    -1.81969   -2.30723   -2.07258   -1.52476    0.20051    0.00867   -1.03723   -0.48754   -0.25289    0.88104
  11    -1.03473   -1.34882   -1.09998   -1.20344    0.22865    0.02686   -0.88936   -0.31408   -0.06525    0.76380

k-point #     3:    0.000000    0.500000    0.500000    0.375000
 state   E_KS       E_HF       E_GW       Sx         Re(Sc)     Im(Sc)     Vxc        DE_HF      DE_GW      Znk
   1   -14.68627  -16.50308  -16.31014   -4.72516    0.17351    0.00002   -2.90835   -1.81681   -1.62387    0.98818
   2   -11.48915  -12.99338  -12.79204   -4.36834    0.18150    0.00003   -2.86410   -1.50424   -1.30290    0.98500
   3   -11.48915  -12.99338  -12.79204   -4.36834    0.18150    0.00003   -2.86410   -1.50424   -1.30289    0.98500
   4   -11.48915  -12.99338  -12.79204   -4.36834    0.18151    0.00003   -2.86410   -1.50424   -1.30289    0.98500
   5    -6.24405   -7.20722   -7.01657   -3.72737    0.17144   -0.00001   -2.76420   -0.96317   -0.77252    0.97574
   6    -6.24403   -7.20721   -7.01657   -3.72738    0.17144   -0.00000   -2.76421   -0.96317   -0.77253    0.97575
   7    -6.24346   -7.20934   -7.01757   -3.73064    0.17268    0.00004   -2.76476   -0.96588   -0.77411    0.97593
   8    -6.24344   -7.20932   -7.01755   -3.73066    0.17268    0.00003   -2.76478   -0.96588   -0.77412    0.97594
   9    -6.24344   -7.20932   -7.01755   -3.73066    0.17268    0.00003   -2.76478   -0.96588   -0.77411    0.97593
  10    -1.81908   -2.30620   -2.07323   -1.52531    0.19899    0.00760   -1.03818   -0.48712   -0.25415    0.88205
  11    -1.03685   -1.35370   -1.10126   -1.20433    0.23184    0.02427   -0.88748   -0.31684   -0.06441    0.75775

  """


def test_k_points_from_evalqp():

    ref = {
        1: {'k_point': [0.000000, 0.000000, 0.000000], 'weight': 0.125000},
        2: {'k_point': [0.000000, 0.000000, 0.500000], 'weight': 0.500000},
        3: {'k_point': [0.000000, 0.500000, 0.500000], 'weight': 0.375000}
        }

    # Note, interaction with the file system should be avoided where possible
    file_name = "EVALQP_TEST.DAT"

    write_file(file_name, evalqp_str)
    k_points_and_weights = k_points_from_evalqp(file_name)
    os.remove(file_name)

    assert k_points_and_weights == ref, "k_points_and_weights should match reference"


def test_n_states_from_evalqp():

    # String defined above
    assert n_states_from_evalqp(evalqp_str) == 11, "Final k-point in string should have 11 states"

    evalqp_str_with_one_kpoint = """k-point #     1:    0.000000    0.000000    0.000000    0.125000
 state   E_KS       E_HF       E_GW       Sx         Re(Sc)     Im(Sc)     Vxc        DE_HF      DE_GW      Znk
   1   -14.68627  -16.50308  -16.31014   -4.72516    0.17351    0.00002   -2.90835   -1.81681   -1.62387    0.98817
   2   -11.48914  -12.99338  -12.79204   -4.36834    0.18150    0.00003   -2.86410   -1.50424   -1.30290    0.98500
   3   -11.48914  -12.99338  -12.79204   -4.36834    0.18150    0.00003   -2.86410   -1.50424   -1.30290    0.98500
   4   -11.48914  -12.99338  -12.79203   -4.36834    0.18151    0.00003   -2.86410   -1.50424   -1.30289    0.98500
   5    -6.24399   -7.20716   -7.01655   -3.72742    0.17144    0.00000   -2.76425   -0.96317   -0.77257    0.97579
   6    -6.24399   -7.20716   -7.01652   -3.72742    0.17144   -0.00001   -2.76425   -0.96317   -0.77253    0.97575
   7    -6.24340   -7.20928   -7.01752   -3.73069    0.17268    0.00003   -2.76481   -0.96588   -0.77412    0.97594
   8    -6.24340   -7.20928   -7.01752   -3.73069    0.17268    0.00003   -2.76481   -0.96588   -0.77412    0.97594
   9    -6.24340   -7.20928   -7.01752   -3.73069    0.17268    0.00003   -2.76481   -0.96588   -0.77412    0.97594

  """

    assert n_states_from_evalqp(evalqp_str_with_one_kpoint) == 9, "k-point 1 should have 9 states"


def test_parse_evalqp():

    # Energies for all states, per k-point. Copied from the reference string above
    energies_1 = np.array([
        [
            -14.68627,
            -16.50308,
            -16.31014,
            -4.72516,
            0.17351,
            0.00002,
            -2.90835,
            -1.81681,
            -1.62387,
            0.98817],
        [
            -11.48914,
            -12.99338,
            -12.79204,
            -4.36834,
            0.18150,
            0.00003,
            -2.86410,
            -1.50424,
            -1.30290,
            0.98500
            ],
        [
            -11.48914,
            -12.99338,
            -12.79204,
            -4.36834,
            0.18150,
            0.00003,
            -2.86410,
            -1.50424,
            -1.30290,
            0.98500
            ],
        [
            -11.48914,
            -12.99338,
            -12.79203,
            -4.36834,
            0.18151,
            0.00003,
            -2.86410,
            -1.50424,
            -1.30289,
            0.98500
            ],
        [
            -6.24399,
            -7.20716,
            -7.01655,
            -3.72742,
            0.17144,
            0.00000,
            -2.76425,
            -0.96317,
            -0.77257,
            0.97579
            ],
        [
            -6.24399,
            -7.20716,
            -7.01652,
            -3.72742,
            0.17144,
            -0.00001,
            -2.76425,
            -0.96317,
            -0.77253,
            0.97575
            ],
        [
            -6.24340,
            -7.20928,
            -7.01752,
            -3.73069,
            0.17268,
            0.00003,
            -2.76481,
            -0.96588,
            -0.77412,
            0.97594
            ],
        [
            -6.24340,
            -7.20928,
            -7.01752,
            -3.73069,
            0.17268,
            0.00003,
            -2.76481,
            -0.96588,
            -0.77412,
            0.97594
            ],
        [
            -6.24340,
            -7.20928,
            -7.01752,
            -3.73069,
            0.17268,
            0.00003,
            -2.76481,
            -0.96588,
            -0.77412,
            0.97594
            ],
        [
            -1.82156,
            -2.31023,
            -2.07248,
            -1.52298,
            0.20180,
            0.00623,
            -1.03431,
            -0.48867,
            -0.25092,
            0.87468
            ],
        [
            -1.00771,
            -1.29282,
            -1.07701,
            -1.21488,
            0.19790,
            0.02221,
            -0.92977,
            -0.28511,
            -0.06930,
            0.79459
            ]])

    energies_2 = np.array([
        [
            -14.68627,
            -16.50308,
            -16.31014,
            -4.72516,
            0.17351,
            0.00002,
            -2.90835,
            -1.81681,
            -1.62387,
            0.98817
        ],
        [
            -11.48915,
            -12.99338,
            -12.79204,
            -4.36834,
            0.18150,
            0.00003,
            -2.86410,
            -1.50424,
            -1.30289,
            0.98500
            ],
        [
            -11.48915,
            -12.99338,
            -12.79204,
            -4.36834,
            0.18151,
            0.00003,
            -2.86410,
            -1.50424,
            -1.30289,
            0.98500
            ],
        [
            -11.48914,
            -12.99338,
            -12.79204,
            -4.36834,
            0.18150,
            0.00003,
            -2.86410,
            -1.50424,
            -1.30290,
            0.98500
            ],
        [
            -6.24401,
            -7.20719,
            -7.01653,
            -3.72740,
            0.17145,
            -0.00001,
            -2.76423,
            -0.96318,
            -0.77252,
            0.97574
            ],
        [
            -6.24401,
            -7.20719,
            -7.01650,
            -3.72740,
            0.17145,
            -0.00001,
            -2.76423,
            -0.96318,
            -0.77249,
            0.97570
            ],
        [
            -6.24345,
            -7.20933,
            -7.01756,
            -3.73065,
            0.17268,
            0.00003,
            -2.76477,
            -0.96588,
            -0.77411,
            0.97593
            ],
        [
            -6.24344,
            -7.20932,
            -7.01754,
            -3.73066,
            0.17268,
            0.00003,
            -2.76478,
            -0.96588,
            -0.77410,
            0.97593
            ],
        [
            -6.24344,
            -7.20932,
            -7.01755,
            -3.73066,
            0.17268,
            0.00003,
            -2.76478,
            -0.96588,
            -0.77411,
            0.97593
            ],
        [
            -1.81969,
            -2.30723,
            -2.07258,
            -1.52476,
            0.20051,
            0.00867,
            -1.03723,
            -0.48754,
            -0.25289,
            0.88104
            ],
        [
            -1.03473,
            -1.34882,
            -1.09998,
            -1.20344,
            0.22865,
            0.02686,
            -0.88936,
            -0.31408,
            -0.06525,
            0.76380
            ]])

    energies_3 = np.array([
        [
            -14.68627,
            -16.50308,
            -16.31014,
            -4.72516,
            0.17351,
            0.00002,
            -2.90835,
            -1.81681,
            -1.62387,
            0.98818
        ],
        [
            -11.48915,
            -12.99338,
            -12.79204,
            -4.36834,
            0.18150,
            0.00003,
            -2.86410,
            -1.50424,
            -1.30290,
            0.98500
            ],
        [
            -11.48915,
            -12.99338,
            -12.79204,
            -4.36834,
            0.18150,
            0.00003,
            -2.86410,
            -1.50424,
            -1.30289,
            0.98500
            ],
        [
            -11.48915,
            -12.99338,
            -12.79204,
            -4.36834,
            0.18151,
            0.00003,
            -2.86410,
            -1.50424,
            -1.30289,
            0.98500
            ],
        [
            -6.24405,
            -7.20722,
            -7.01657,
            -3.72737,
            0.17144,
            -0.00001,
            -2.76420,
            -0.96317,
            -0.77252,
            0.97574
            ],
        [
            -6.24403,
            -7.20721,
            -7.01657,
            -3.72738,
            0.17144,
            -0.00000,
            -2.76421,
            -0.96317,
            -0.77253,
            0.97575
            ],
        [
            -6.24346,
            -7.20934,
            -7.01757,
            -3.73064,
            0.17268,
            0.00004,
            -2.76476,
            -0.96588,
            -0.77411,
            0.97593
            ],
        [
            -6.24344,
            -7.20932,
            -7.01755,
            -3.73066,
            0.17268,
            0.00003,
            -2.76478,
            -0.96588,
            -0.77412,
            0.97594
            ],
        [
            -6.24344,
            -7.20932,
            -7.01755,
            -3.73066,
            0.17268,
            0.00003,
            -2.76478,
            -0.96588,
            -0.77411,
            0.97593
            ],
        [
            -1.81908,
            -2.30620,
            -2.07323,
            -1.52531,
            0.19899,
            0.00760,
            -1.03818,
            -0.48712,
            -0.25415,
            0.88205
            ],
        [
            -1.03685,
            -1.35370,
            -1.10126,
            -1.20433,
            0.23184,
            0.02427,
            -0.88748,
            -0.31684,
            -0.06441,
            0.75775
            ]])

    ref = {
        1: {'k_point': [0.000000, 0.000000, 0.000000], 'weight': 0.125000},
        2: {'k_point': [0.000000, 0.000000, 0.500000], 'weight': 0.500000},
        3: {'k_point': [0.000000, 0.500000, 0.500000], 'weight': 0.375000}
        }

    # Note, interaction with the file system should be avoided where possible
    file_name = "EVALQP_TEST.DAT"
    write_file(file_name, evalqp_str)
    output = parse_evalqp(file_name)
    os.remove(file_name)

    # k-points
    assert len(output) == 3, "Expect 3 k-points"
    assert min([ik for ik in output.keys()]) == 1, 'k-point indexing starts at 1'
    assert output[1]['k_point'] == ref[1]['k_point'], "Compare k-point 1 to reference"
    assert output[2]['k_point'] == ref[2]['k_point'], "Compare k-point 2 to reference"
    assert output[3]['k_point'] == ref[3]['k_point'], "Compare k-point 3 to reference"

    # Weights
    assert output[1]['weight'] == ref[1]['weight'], "Compare weight 1 to reference"
    assert output[2]['weight'] == ref[2]['weight'], "Compare weight 2 to reference"
    assert output[3]['weight'] == ref[3]['weight'], "Compare weight 3 to reference"

    # Energies
    assert output[1]['energies'].shape == (11, 10), 'rows = 11 states and cols = 10 energies'
    assert np.allclose(output[1]['energies'], energies_1), "Compare energies 1 to reference"
    assert np.allclose(output[2]['energies'], energies_2), "Compare energies 2 to reference"
    assert np.allclose(output[3]['energies'], energies_3), "Compare energies 3 to reference"


# Example output from VXCNN.DAT
vxc_string = """ik=   1    vkl=  0.0000  0.0000  0.0000
   1       -2.908349       -0.000000
   2       -2.864103        0.000000
   3       -2.864103       -0.000000
   4       -2.864103       -0.000000
   5       -2.764246       -0.000000
   6       -2.764246        0.000000
   7       -2.764809       -0.000000
   8       -2.764809       -0.000000
   9       -2.764809       -0.000000
  10       -1.034312       -0.000000
  11       -0.929773       -0.000000

ik=   2    vkl=  0.0000  0.0000  0.5000
 1       -2.908349        0.000000
 2       -2.864100       -0.000000
 3       -2.864100        0.000000
 4       -2.864101       -0.000000
 5       -2.764227       -0.000000
 6       -2.764227        0.000000
 7       -2.764770        0.000000
 8       -2.764777        0.000000
 9       -2.764777        0.000000
10       -1.037228       -0.000000
11       -0.889360        0.000000

ik=   3    vkl=  0.0000  0.5000  0.5000
   1       -2.908349       -0.000000
   2       -2.864099       -0.000000
   3       -2.864099       -0.000000
   4       -2.864099        0.000000
   5       -2.764195        0.000000
   6       -2.764208        0.000000
   7       -2.764760       -0.000000
   8       -2.764780       -0.000000
   9       -2.764780       -0.000000
  10       -1.038185        0.000000
  11       -0.887485       -0.000000
"""


def test_vkl_from_vxc():

    # Note, interaction with the file system should be avoided where possible
    file_name = "VXCNN_TEST.DAT"
    vkl_ref = {
        1: [0.0000, 0.0000, 0.0000], 2: [0.0000, 0.0000, 0.5000], 3: [0.0000, 0.5000, 0.5000]
        }

    write_file(file_name, vxc_string)
    output = vkl_from_vxc(file_name)
    os.remove(file_name)

    assert len(output) == 3, "Expect 3 vkl/k-points"
    assert output == vkl_ref, "vkl values equal to vkl_ref"


def test_parse_vxcnn():

    # Reference V_xc extracted from vxc_string, defined above
    v_xc_1 = np.array([[-2.908349, -0.000000], [-2.864103, 0.000000], [-2.864103, -0.000000],
                       [-2.864103, -0.000000], [-2.764246, -0.000000], [-2.764246, 0.000000],
                       [-2.764809, -0.000000], [-2.764809, -0.000000], [-2.764809, -0.000000],
                       [-1.034312, -0.000000], [-0.929773, -0.000000]])

    v_xc_2 = np.array([[-2.908349, 0.000000], [-2.864100, -0.000000], [-2.864100, 0.000000],
                       [-2.864101, -0.000000], [-2.764227, -0.000000], [-2.764227, 0.000000],
                       [-2.764770, 0.000000], [-2.764777, 0.000000], [-2.764777, 0.000000],
                       [-1.037228, -0.000000], [-0.889360, 0.000000]])

    v_xc_3 = np.array([[-2.908349, -0.000000], [-2.864099, -0.000000], [-2.864099, -0.000000],
                       [-2.864099, 0.000000], [-2.764195, 0.000000], [-2.764208, 0.000000],
                       [-2.764760, -0.000000], [-2.764780, -0.000000], [-2.764780, -0.000000],
                       [-1.038185, 0.000000], [-0.887485, -0.000000]])

    # Note, interaction with the file system should be avoided where possible
    file_name = "VXCNN_TEST.DAT"
    write_file(file_name, vxc_string)
    output = parse_vxcnn(file_name)
    os.remove(file_name)

    assert [key for key in output[1].keys()] == ['vkl', 'v_xc_nn'], "Key consistency for ik=1 of parsed vxcnn"
    assert [key for key in output[2].keys()] == ['vkl', 'v_xc_nn'], "Key consistency for ik=2 of parsed vxcnn"
    assert [key for key in output[3].keys()] == ['vkl', 'v_xc_nn'], "Key consistency for ik=3 of parsed vxcnn"

    assert output[1]['vkl'] == [0.0000, 0.0000, 0.0000], "vkl (ik=1)"
    assert output[2]['vkl'] == [0.0000, 0.0000, 0.5000], "vkl (ik=2)"
    assert output[3]['vkl'] == [0.0000, 0.5000, 0.5000], "vkl (ik=3)"

    assert output[1]['v_xc_nn'].shape == (11, 2), "Expect V_xc to have 2 cols for 11 states"
    assert output[2]['v_xc_nn'].shape == (11, 2), "Expect V_xc to have 2 cols for 11 states"
    assert output[3]['v_xc_nn'].shape == (11, 2), "Expect V_xc to have 2 cols for 11 states"

    assert np.allclose(output[1]['v_xc_nn'], v_xc_1), "v_xc_nn for ik=1"
    assert np.allclose(output[2]['v_xc_nn'], v_xc_2), "v_xc_nn for ik=2"
    assert np.allclose(output[3]['v_xc_nn'], v_xc_3), "v_xc_nn for ik=3"


# Example output from EPS00_GW.OUT
# File bizarrely begins with a blank line
eps_string = """
 (dielectric tensor, random phase approximation)

 frequency index and value:      1    0.00529953
 real part, imaginary part below
    8.31881773    0.00000000    0.00000000         0.00000000    0.00000000    0.00000000
    0.00000000    8.31881773    0.00000000         0.00000000    0.00000000    0.00000000
    0.00000000    0.00000000    8.31881773         0.00000000    0.00000000    0.00000000

 frequency index and value:      2    0.02771249
 real part, imaginary part below
    8.22189228    0.00000000    0.00000000        -0.00000000    0.00000000    0.00000000
    0.00000000    8.22189228    0.00000000         0.00000000   -0.00000000    0.00000000
    0.00000000    0.00000000    8.22189228         0.00000000    0.00000000   -0.00000000


    """


def test_parse_eps00_frequencies():

    eps_string2 = """
 (dielectric tensor, random phase approximation)

 frequency index and value:      1    0.00529953
 real part, imaginary part below
    8.31881773    0.00000000    0.00000000         0.00000000    0.00000000    0.00000000
    0.00000000    8.31881773    0.00000000         0.00000000    0.00000000    0.00000000
    0.00000000    0.00000000    8.31881773         0.00000000    0.00000000    0.00000000

 frequency index and value:      2    0.02771249
 real part, imaginary part below
    8.22189228    0.00000000    0.00000000        -0.00000000    0.00000000    0.00000000
    0.00000000    8.22189228    0.00000000         0.00000000   -0.00000000    0.00000000
    0.00000000    0.00000000    8.22189228         0.00000000    0.00000000   -0.00000000

 frequency index and value:      3    0.06718440
 real part, imaginary part below
    7.78004308    0.00000000    0.00000000         0.00000000    0.00000000    0.00000000
    0.00000000    7.78004308    0.00000000         0.00000000    0.00000000    0.00000000
    0.00000000    0.00000000    7.78004308         0.00000000    0.00000000    0.00000000

    """

    line = get_new_line_indices(eps_string2)
    assert eps_string2[line[0]:line[1]].isspace(), "First line of eps_string2 must be a whiteline"

    ref = {1: 0.00529953, 2: 0.02771249, 3: 0.06718440}
    assert parse_eps00_frequencies(eps_string2) == ref, "Frequency grid for eps00"


def test_parse_eps00_gw():

    line = get_new_line_indices(eps_string)
    assert eps_string[line[0]:line[1]].isspace(), "First line of eps_string must be a whiteline"

    ref = {
        1: {
            'frequency': 0.00529953,
            'eps00': {
                're': np.array([[8.31881773, 0., 0.], [0., 8.31881773, 0.], [0., 0., 8.31881773]]),
                'img': np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
                }
            },
        2: {
            'frequency': 0.02771249,
            'eps00': {
                're': np.array([[8.22189228, 0., 0.], [0., 8.22189228, 0.], [0., 0., 8.22189228]]),
                'img': np.array([[-0., 0., 0.], [0., -0., 0.], [0., 0., -0.]])
                }
            }
        }

    output = parse_eps00_gw(eps_string)

    assert len(output) == 2, "2 frequency points"
    assert [k for k in output[1].keys()] == ['frequency', 'eps00'], "Frequency point 1 keys "
    assert [k for k in output[2].keys()] == ['frequency', 'eps00'], "Frequency point 2 keys "

    assert output[1]['frequency'] == 0.00529953, "Frequency point 1 value"
    assert output[2]['frequency'] == 0.02771249, "Frequency point 2 value"

    assert np.allclose(output[1]['eps00']['re'], ref[1]['eps00']['re']),"Re{eps00} at frequency point 1"
    assert np.allclose(output[1]['eps00']['img'], ref[1]['eps00']['img']),"Im{eps00} at frequency point 1"

    assert np.allclose(output[2]['eps00']['re'], ref[2]['eps00']['re']),"Re{eps00} at frequency point 2"
    assert np.allclose(output[2]['eps00']['img'], ref[2]['eps00']['img']),"Im{eps00} at frequency point 2"
