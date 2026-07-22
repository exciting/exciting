import os
from typing import Tuple

import excitingscripts
import numpy as np
import pytest
from excitingscripts.execute.convergence_test import execute_convergence_test
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def info_out_mock(tmp_path) -> Tuple[MockFile, MockFile, MockFile, MockFile]:
    """ Mock 'INFO.OUT' data.
    """
    info_out_4_5_str = """
    ================================================================================
    | EXCITING OXYGEN started                                                      =
    | version hash id: e93826b67f8bc6368a11d0f93b8a28d9cbd27476                    =
    |                                                                              =
    | compiler: ifort (IFORT) 2021.6.0 20220226                                    =
    |                                                                              =
    |                                                                              =
    | Date (DD-MM-YYYY) : 07-12-2022                                               =
    | Time (hh:mm:ss)   : 14:16:33                                                 =
    |                                                                              =
    | All units are atomic (Hartree, Bohr, etc.)                                   =
    ================================================================================
     
    ********************************************************************************
    * Ground-state run starting from atomic densities                              *
    ********************************************************************************
     
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    + Starting initialization                                                      +
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     
     Lattice vectors (cartesian) :
          3.8600500000      3.8600500000      0.0000000000
          3.8600500000      0.0000000000      3.8600500000
          0.0000000000      3.8600500000      3.8600500000
     
     Reciprocal lattice vectors (cartesian) :
          0.8138735647      0.8138735647     -0.8138735647
          0.8138735647     -0.8138735647      0.8138735647
         -0.8138735647      0.8138735647      0.8138735647
     
     Unit cell volume                           :     115.0293819379
     Brillouin zone volume                      :       2.1564074262
     
     Species :    1 (Ag)
         parameters loaded from                 :    Ag.xml
         name                                   :    silver
         nuclear charge                         :     -47.00000000
         electronic charge                      :      47.00000000
         atomic mass                            :  196631.69970000
         muffin-tin radius                      :       2.00000000
         # of radial points in muffin-tin       :     400
     
         atomic positions (lattice) :
           1 :   0.00000000  0.00000000  0.00000000
     
     Total number of atoms per unit cell        :       1
     
     Spin treatment                             :    spin-unpolarised
     
     Number of Bravais lattice symmetries       :      48
     Number of crystal symmetries               :      48
     
     k-point grid                               :       4    4    4
     Total number of k-points                   :       8
     k-point set is reduced with crystal symmetries
     
     R^MT_min * |G+k|_max (rgkmax)              :       5.00000000
     Species with R^MT_min                      :       1 (Ag)
     Maximum |G+k| for APW functions            :       2.50000000
     Maximum |G| for potential and density      :      12.00000000
     
     G-vector grid sizes                        :      24    24    24
     Total number of G-vectors                  :    3383
     
     Maximum angular momentum used for
         APW functions                          :       8
         computing H and O matrix elements      :       8
         potential and density                  :       8
         inner part of muffin-tin               :       2
     
     Total nuclear charge                       :     -47.00000000
     Total electronic charge                    :      47.00000000
     Total core charge                          :      28.00000000
     Total valence charge                       :      19.00000000
     
     Number of empty states                     :       5
     Total number of valence states             :      15
     
     Maximum Hamiltonian size                   :      47
     Maximum number of plane-waves              :      34
     Total number of local-orbitals             :      13
     
     Exchange-correlation type                  :      22
         PBEsol, Phys. Rev. Lett. 100, 136406 (2008)
         Generalised gradient approximation (GGA)
     
     Smearing scheme                            :    Gaussian
     Smearing width                             :       0.00100000
     
     Using multisecant Broyden potential mixing
     
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    + Ending initialization                                                        +
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     
    ********************************************************************************
    * Groundstate module started                                                   *
    ********************************************************************************
     Output level for this task is set to normal
     
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    + Self-consistent loop started                                                 +
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     Density and potential initialised from atomic data
     
     
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    + SCF iteration number :   14                                                  +
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     Total energy                               :     -5314.29682116
     _______________________________________________________________
     Fermi energy                               :         0.33193980
     Kinetic energy                             :      5561.19075068
     Coulomb energy                             :    -10728.57982545
     Exchange energy                            :      -144.05330541
     Correlation energy                         :        -2.85444098
     
     DOS at Fermi energy (states/Ha/cell)       :        94.47871339
     
     Electron charges :
         core                                   :        28.00000000
         core leakage                           :         0.00000062
         valence                                :        19.00000000
         interstitial                           :         3.02377384
         charge in muffin-tin spheres :
                      atom     1    Ag          :        43.97622616
         total charge in muffin-tins            :        43.97622616
         total charge                           :        47.00000000
     
     Wall time (seconds)                        :        58.60
     
     RMS change in effective potential (target) :  0.597519E-10  ( 0.100000E-05)
     Absolute change in total energy   (target) :  0.305927E-07  ( 0.100000E-05)
     Charge distance                   (target) :  0.349170E-10  ( 0.100000E-04)
                                                                                    
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    | Convergency criteria checked for the last 2 iterations                       +
    | Convergence targets achieved. Performing final SCF iteration                 +
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     Total energy                               :     -5314.29682116
     _______________________________________________________________
     Fermi energy                               :         0.33193980
     Kinetic energy                             :      5561.19075068
     Coulomb energy                             :    -10728.57982545
     Exchange energy                            :      -144.05330541
     Correlation energy                         :        -2.85444098
     
     DOS at Fermi energy (states/Ha/cell)       :        94.47871320
     
     Electron charges :
         core                                   :        28.00000000
         core leakage                           :         0.00000062
         valence                                :        19.00000000
         interstitial                           :         3.02377384
         charge in muffin-tin spheres :
                      atom     1    Ag          :        43.97622616
         total charge in muffin-tins            :        43.97622616
         total charge                           :        47.00000000
     
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    + Self-consistent loop stopped                                                 +
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     STATE.OUT is written
     
    ********************************************************************************
    * Groundstate module stopped                                                   *
    ********************************************************************************
     
     Total time spent (seconds)                 :        51.62
    ================================================================================
    | EXCITING OXYGEN stopped                                                      =
    ================================================================================
    """
    info_out_4_6_str = info_out_4_5_str.replace('-5314.29682116', '-5314.65194126')
    info_out_6_5_str = info_out_4_5_str.replace('-5314.29682116', '-5314.3178557')
    info_out_6_6_str = info_out_4_5_str.replace('-5314.29682116', '-5314.65182124')

    (k_initial, k_final, dk) = (4, 7, 2)
    rgkmax_initial = 5
    rgkmax_final = 7
    for k in range(k_initial, k_final, dk):
        for rgkmax in range(rgkmax_initial, rgkmax_final):
            os.makedirs(os.path.dirname(tmp_path / f"{k}_{rgkmax}/INFO.OUT"), exist_ok=True)

    info_out_4_5_file = tmp_path / "4_5/INFO.OUT"
    info_out_4_5_file.write_text(info_out_4_5_str)

    info_out_4_6_file = tmp_path / "4_6/INFO.OUT"
    info_out_4_6_file.write_text(info_out_4_6_str)

    info_out_6_5_file = tmp_path / "6_5/INFO.OUT"
    info_out_6_5_file.write_text(info_out_6_5_str)

    info_out_6_6_file = tmp_path / "6_6/INFO.OUT"
    info_out_6_6_file.write_text(info_out_6_6_str)

    return MockFile(info_out_4_5_file, info_out_4_5_str), MockFile(info_out_4_6_file, info_out_4_6_str), \
           MockFile(info_out_6_5_file, info_out_6_5_str), MockFile(info_out_6_6_file, info_out_6_6_str)


def test_execute_convergence_test(monkeypatch, info_out_mock, tmp_path):
    number_calculations = 4

    def replace_run_exciting(*args, **kwargs):
        for i in range(number_calculations):
            with open(info_out_mock[i].full_path, "w") as f:
                f.write(info_out_mock[i].string)

    monkeypatch.setattr(excitingscripts.execute.convergence_test, 'run_exciting', replace_run_exciting)

    convergence_test_ref =  np.array([[4, 5, -5314.29682116],
                                      [4, 6, -5314.65194126],
                                      [6, 5, -5314.3178557],
                                      [6, 6, -5314.65182124]])

    assert np.allclose(execute_convergence_test(4, 6, 2, 5, 6, 1, root_directory=tmp_path), convergence_test_ref)

