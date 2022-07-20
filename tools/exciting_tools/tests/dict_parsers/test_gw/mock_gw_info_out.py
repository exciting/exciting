"""Outputs from GW_INFO.OUT
"""
import pytest
from excitingtools.utils.test_utils import MockFile

from excitingtools.exciting_dict_parsers.gw_info_parser import _file_name


@pytest.fixture
def zro2_gw_info_out_mock(tmp_path):
    file = tmp_path / _file_name
    file.write_text(zro2_gw_info_out)
    return MockFile(file, zro2_gw_info_out)


@pytest.fixture
def si_2_gw_info_out_mock(tmp_path):
    file = tmp_path / _file_name
    file.write_text(si_2_gw_info_out)
    return MockFile(file, si_2_gw_info_out)


# GW output for ZrO2. This reports an indirect gap, as well as its direct gap
zro2_gw_info_out = """
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
   imfreq : weights calculated for imaginary frequencies

--------------------------------------------------------------------------------

 Correlation self-energy parameters:
 Solution of the QP equation:
   0 - perturbative solution
 Energy alignment:
   0 - no alignment
 Analytic continuation method:
  PADE - Thiele's reciprocal difference method
  Citation: H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179 (1977)
 Scheme to treat singularities:
  Auxiliary function method "mpb"
  Citation: S. Massidda, M. Posternak, and A. Baldereschi, PRB 48, 5058 (1993)

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

# Reports a direct gap only
si_2_gw_info_out = """
================================================================================
=                            GW input parameters                               =
================================================================================


 GW taskname:

   g0w0 - G0W0 run

--------------------------------------------------------------------------------

 Frequency integration parameters:
 Number of frequencies:           16
 Cutoff frequency:    1.0000000000000000     
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
     Angular momentum cutoff:            3
     Linear dependence tolerance factor:    1.0000000000000000E-004
   Interstitial:
     Plane wave cutoff (in units of Gkmax):    1.0000000000000000     

--------------------------------------------------------------------------------

 Bare Coulomb potential parameters:
   Plane wave cutoff (in units of Gkmax*input%gw%MixBasis%gmb):    2.0000000000000000     
   Error tolerance for structure constants:    1.0000000000000001E-015
   Tolerance factor to reduce the MB size based on
   the eigenvectors of the bare Coulomb potential:   0.10000000000000001     

--------------------------------------------------------------------------------

 Screened Coulomb potential:
   Full-frequency Random-Phase Approximation

--------------------------------------------------------------------------------

 Core electrons treatment:
   all - Core states are included in all calculations

--------------------------------------------------------------------------------

 Interval of quasiparticle states (ibgw, nbgw):       1     15

 Number of empty states (GW):          100

 k/q-points grid:            1           1           1

--------------------------------------------------------------------------------
-                           Mixed product WF info                              -
--------------------------------------------------------------------------------

  Maximal number of MT wavefunctions per atom:          159
  Total number of MT wavefunctions:                     318
  Maximal number of PW wavefunctions:                   169
  Total number of mixed-product wavefunctions:          487


--------------------------------------------------------------------------------
-                               frequency grid                                 -
--------------------------------------------------------------------------------

 Type: < fgrid > gauleg2                                 
 Frequency axis: < fconv > imfreq                                  
 Number of frequencies: < nomeg >           16
 Cutoff frequency: < freqmax >    1.0000000000000000     
 frequency list: < #    freqs    weight > 
   1  1.9855071751E-02  5.0614268145E-02
   2  0.1016667613      0.1111905172    
   3  0.2372337950      0.1568533229    
   4  0.4082826788      0.1813418917    
   5  0.5917173212      0.1813418917    
   6  0.7627662050      0.1568533229    
   7  0.8983332387      0.1111905172    
   8  0.9801449282      5.0614268145E-02
   9   1.020257282      5.2685653046E-02
  10   1.113172659      0.1377821040    
  11   1.311017706      0.2695943819    
  12   1.689996159      0.5179282224    
  13   2.449283430       1.087868072    
  14   4.215251035       2.787023374    
  15   9.836056419       10.75746081    
  16   50.36496531       128.3896574    

--------------------------------------------------------------------------------
-                       Kohn-Sham eigenstates summary                          -
--------------------------------------------------------------------------------

 Maximum number of LAPW states:                      177
 Minimal number of LAPW states:                      177
 Number of states used in GW:
     - total KS                                      105
     - occupied                                        4
     - unoccupied                                    100
     - dielectric function                           105
     - self energy                                   105
 Energy of the highest unoccupied state:        3.407168
 Number of valence electrons:                          8
 Number of valence electrons treated in GW:            8

--------------------------------------------------------------------------------
-                          Kohn-Sham band structure                            -
--------------------------------------------------------------------------------

 Fermi energy:     0.0000
 Energy range:    -0.4868    3.4072
 Band index of VBM:   4
 Band index of CBm:   5

 Direct BandGap (eV):                      2.3903
 at k      =    0.000   0.000   0.000 ik =     1

================================================================================
=                                  GW cycle                                    =
================================================================================

 (task_gw): q-point cycle, iq =            1

--------------------------------------------------------------------------------
-                            G0W0 band structure                               -
--------------------------------------------------------------------------------

 Fermi energy:     0.0176
 Energy range:    -0.4799    0.5045
 Band index of VBM:   4
 Band index of CBm:   5

 Direct BandGap (eV):                      3.2457
 at k      =    0.000   0.000   0.000 ik =     1

================================================================================
=                          GW timing info (seconds)                            =
================================================================================

 Initialization                             :         0.44
     - init_scf                             :         0.31
     - init_kpt                             :         0.01
     - init_eval                            :         0.00
     - init_freq                            :         0.00
     - init_mb                              :         0.05
 Subroutines                                : 
     - calcpmat                             :         0.39
     - calcbarcmb                           :         0.29
     - BZ integration weights               :         0.02
     Dielectric function                    :         0.55
     - head                                 :         0.00
     - wings                                :         0.00
     - body                                 :         0.00
     - inversion                            :         0.27
     WF products expansion                  :         0.03
     - diagsgi                              :         0.02
     - calcmpwipw                           :         0.01
     - calcmicm                             :         0.00
     - calcminc                             :         0.00
     - calcminm                             :         0.00
     Self-energy                            :         2.04
     - calcselfx                            :         0.04
     - calcselfc                            :         2.00
     - calcvxcnn                            :         0.11
     - input/output                         :         0.00
 _________________________________________________________
 Total                                      :         4.25
"""

