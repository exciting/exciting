import os
from typing import Tuple

import excitingscripts
import numpy as np
import pytest
from excitingscripts.execute.volume_optimization import execute_volume_optimization
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def input_xml_mock(tmp_path) -> MockFile:
    """ Mock 'input.xml' data.
    """
    input_xml_str = """<?xml version="1.0" encoding="UTF-8"?>
    <input>

       <title>Silver</title>
    
       <structure speciespath="$EXCITINGROOT/species">
    
          <crystal scale="7.729">
             <basevect> 0.5     0.5     0.0 </basevect>
             <basevect> 0.5     0.0     0.5 </basevect>
             <basevect> 0.0     0.5     0.5 </basevect>
          </crystal>
    
          <species speciesfile="Ag.xml">
             <atom coord="0.00 0.00 0.00" />
          </species>
    
       </structure>
    
       <groundstate
          ngridk="8 8 8"
          swidth="0.01"
          rgkmax="7.5"
          xctype="GGA_PBE_SOL">
       </groundstate>

    </input>
    """

    input_xml_file = tmp_path / "input.xml"
    input_xml_file.write_text(input_xml_str)

    return MockFile(input_xml_file, input_xml_str)


@pytest.fixture
def info_out_mock(tmp_path) -> Tuple[MockFile, MockFile, MockFile]:
    """ Mock 'INFO.OUT' data.
    """
    info_out_vol_1_str = """
    ================================================================================
    | EXCITING OXYGEN started                                                      =
    | version hash id: e93826b67f8bc6368a11d0f93b8a28d9cbd27476                    =
    |                                                                              =
    | compiler: ifort (IFORT) 2021.6.0 20220226                                    =
    |                                                                              =
    |                                                                              =
    | Date (DD-MM-YYYY) : 15-01-2023                                               =
    | Time (hh:mm:ss)   : 20:59:32                                                 =
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
          3.6712750000      3.6712750000      0.0000000000
          3.6712750000      0.0000000000      3.6712750000
          0.0000000000      3.6712750000      3.6712750000
     
     Reciprocal lattice vectors (cartesian) :
          0.8557225088      0.8557225088     -0.8557225088
          0.8557225088     -0.8557225088      0.8557225088
         -0.8557225088      0.8557225088      0.8557225088
     
     Unit cell volume                           :      98.9647988854
     Brillouin zone volume                      :       2.5064489216
     
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
     
     k-point grid                               :       8    8    8
     Total number of k-points                   :      29
     k-point set is reduced with crystal symmetries
     
     R^MT_min * |G+k|_max (rgkmax)              :       7.50000000
     Species with R^MT_min                      :       1 (Ag)
     Maximum |G+k| for APW functions            :       3.75000000
     Maximum |G| for potential and density      :      12.00000000
     
     G-vector grid sizes                        :      20    20    20
     Total number of G-vectors                  :    2891
     
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
     
     Maximum Hamiltonian size                   :     110
     Maximum number of plane-waves              :      97
     Total number of local-orbitals             :      13
     
     Exchange-correlation type                  :      22
         PBEsol, Phys. Rev. Lett. 100, 136406 (2008)
         Generalised gradient approximation (GGA)
     
     Smearing scheme                            :    Gaussian
     Smearing width                             :       0.01000000
     
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
     Total energy                               :     -5314.77888967
     _______________________________________________________________
     Fermi energy                               :         0.42997785
     Kinetic energy                             :      5560.34911837
     Coulomb energy                             :    -10728.14761006
     Exchange energy                            :      -144.11965057
     Correlation energy                         :        -2.86074741
     
     DOS at Fermi energy (states/Ha/cell)       :        12.55965962
     
     Electron charges :
         core                                   :        28.00000000
         core leakage                           :         0.00000062
         valence                                :        19.00000000
         interstitial                           :         2.65685214
         charge in muffin-tin spheres :
                      atom     1    Ag          :        44.34314786
         total charge in muffin-tins            :        44.34314786
         total charge                           :        47.00000000
     
     Wall time (seconds)                        :        55.95
     
     RMS change in effective potential (target) :  0.902238E-11  ( 0.100000E-05)
     Absolute change in total energy   (target) :  0.438467E-08  ( 0.100000E-05)
     Charge distance                   (target) :  0.312253E-11  ( 0.100000E-04)
                                                                                    
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    | Convergency criteria checked for the last 2 iterations                       +
    | Convergence targets achieved. Performing final SCF iteration                 +
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     Total energy                               :     -5314.77888967
     _______________________________________________________________
     Fermi energy                               :         0.42997785
     Kinetic energy                             :      5560.34911837
     Coulomb energy                             :    -10728.14761006
     Exchange energy                            :      -144.11965057
     Correlation energy                         :        -2.86074741
     
     DOS at Fermi energy (states/Ha/cell)       :        12.55965962
     
     Electron charges :
         core                                   :        28.00000000
         core leakage                           :         0.00000062
         valence                                :        19.00000000
         interstitial                           :         2.65685214
         charge in muffin-tin spheres :
                      atom     1    Ag          :        44.34314786
         total charge in muffin-tins            :        44.34314786
         total charge                           :        47.00000000
     
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    + Self-consistent loop stopped                                                 +
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     STATE.OUT is written
     
    ********************************************************************************
    * Groundstate module stopped                                                   *
    ********************************************************************************
     
     Total time spent (seconds)                 :        48.56
    ================================================================================
    | EXCITING OXYGEN stopped                                                      =
    ================================================================================
    """
    info_out_vol_2_str = info_out_vol_1_str.replace('-5314.77888967', '-5314.78069599') 
    info_out_vol_2_str = info_out_vol_1_str.replace('98.9647988854', '115.42767037')
    info_out_vol_3_str = info_out_vol_1_str.replace('-5314.77888967', '-5314.78200255')
    info_out_vol_3_str = info_out_vol_1_str.replace('98.9647988854', '133.62195691')

    nr_vol = 3
    for i_v in range(0, nr_vol):
        os.makedirs(os.path.dirname(tmp_path / f"volume-{i_v + 1}/INFO.OUT"), exist_ok=True)

    info_out_vol_1_file = tmp_path / "volume-1/INFO.OUT"
    info_out_vol_1_file.write_text(info_out_vol_1_str)

    info_out_vol_2_file = tmp_path / "volume-2/INFO.OUT"
    info_out_vol_2_file.write_text(info_out_vol_2_str)

    info_out_vol_3_file = tmp_path / "volume-3/INFO.OUT"
    info_out_vol_3_file.write_text(info_out_vol_3_str)

    return MockFile(info_out_vol_1_file, info_out_vol_1_str), MockFile(info_out_vol_2_file, info_out_vol_2_str), \
           MockFile(info_out_vol_3_file, info_out_vol_3_str)


def test_execute_volume_optimization(monkeypatch, info_out_mock, input_xml_mock, tmp_path):
    number_calculations = 3

    def replace_run_exciting(*args, **kwargs):
        for i in range(number_calculations):
            with open(info_out_mock[i].full_path, "w") as f:
                f.write(info_out_mock[i].string)

    monkeypatch.setattr(excitingscripts.execute.volume_optimization, 'run_exciting', replace_run_exciting)

    volume_optimization_ref = np.array([[98.96479889, -5314.77888967],
                                        [115.42767037, -5314.78069599],
                                        [133.62195691, -5314.78200255]])

    assert np.allclose(execute_volume_optimization(input_xml_mock.full_path, 3, root_directory=tmp_path),
                       volume_optimization_ref)

