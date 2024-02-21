import os
from typing import Tuple

import excitingscripts
import numpy as np
import pytest
from excitingscripts.execute.elastic_strain import execute_elastic_strain
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def info_out_mock(tmp_path) -> Tuple[MockFile, MockFile, MockFile, MockFile, MockFile]:
    """ Mock 'INFO.OUT' data.
    """
    info_out_strain_1_str = """
    ================================================================================
    | EXCITING NEON started                                                        =
    | version hash id: 5d3aaf6117e024ca12df2ca191030fdb29045cc2                    =
    |                                                                              =
    | compiler: ifort (IFORT) 2021.3.0 20210609                                    =
    |                                                                              =
    |                                                                              =
    | Date (DD-MM-YYYY) : 17-10-2023                                               =
    | Time (hh:mm:ss)   : 16:38:05                                                 =
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
          2.3243679000      4.0259232984      0.0000000000
          4.6487358000      0.0000000000      0.0000000000
          0.0000000000      0.0000000000     10.0000000000
     
     Reciprocal lattice vectors (cartesian) :
          0.0000000000      1.5606818217      0.0000000000
          1.3515901048     -0.7803409108      0.0000000000
          0.0000000000      0.0000000000      0.6283185307
     
     Unit cell volume                           :     187.1545376511
     Brillouin zone volume                      :       1.3253764325
     
     Species :    1 (C)
         parameters loaded from                 :    C.xml
         name                                   :    carbon
         nuclear charge                         :      -6.00000000
         electronic charge                      :       6.00000000
         atomic mass                            :   21894.16673000
         muffin-tin radius                      :       1.20000000
         # of radial points in muffin-tin       :     250
     
         atomic positions (lattice) :
           1 :   0.00000000  0.00000000  0.00000000
           2 :   0.00000000  0.00000000  0.50000000
           3 :   0.66666667  0.66666667  0.00000000
           4 :   0.33333333  0.33333333  0.50000000
     
     Total number of atoms per unit cell        :       4
     
     Spin treatment                             :    spin-unpolarised
     
     Number of Bravais lattice symmetries       :      24
     Number of crystal symmetries               :      24
     
     k-point grid                               :      10   10    4
     Total number of k-points                   :      42
     k-point set is reduced with crystal symmetries
     
     R^MT_min * |G+k|_max (rgkmax)              :       6.00000000
     Species with R^MT_min                      :       1 (C)
     Maximum |G+k| for APW functions            :       5.00000000
     Maximum |G| for potential and density      :      20.00000000
     
     G-vector grid sizes                        :      30    30    64
     Total number of G-vectors                  :   25281
     
     Maximum angular momentum used for
         APW functions                          :       8
         computing H and O matrix elements      :       8
         potential and density                  :       8
         inner part of muffin-tin               :       2
     
     Total nuclear charge                       :     -24.00000000
     Total electronic charge                    :      24.00000000
     Total core charge                          :       8.00000000
     Total valence charge                       :      16.00000000
     
     Number of empty states                     :       5
     Total number of valence states             :      14
     
     Maximum Hamiltonian size                   :     428
     Maximum number of plane-waves              :     412
     Total number of local-orbitals             :      16
     
     Exchange-correlation type                  :      20
         Perdew-Burke-Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)
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
    + SCF iteration number :   13                                                  +
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     Total energy                               :      -152.38457373
     _______________________________________________________________
     Fermi energy                               :         0.33318841
     Kinetic energy                             :       152.06135670
     Coulomb energy                             :      -282.76363837
     Exchange energy                            :       -20.81913392
     Correlation energy                         :        -0.86315814
     
     DOS at Fermi energy (states/Ha/cell)       :        32.73559610
     
     Electron charges :
         core                                   :         8.00000000
         core leakage                           :         0.00261990
         valence                                :        16.00000000
         interstitial                           :         9.41909972
         charge in muffin-tin spheres :
                      atom     1     C          :         3.64065647
                      atom     2     C          :         3.64065647
                      atom     3     C          :         3.64979367
                      atom     4     C          :         3.64979367
         total charge in muffin-tins            :        14.58090028
         total charge                           :        24.00000000
     
     Wall time (seconds)                        :        48.23
     
     RMS change in effective potential (target) :  0.137097E-07  ( 0.100000E-05)
     Absolute change in total energy   (target) :  0.109497E-08  ( 0.100000E-05)
     Charge distance                   (target) :  0.780373E-08  ( 0.100000E-04)
                                                                                    
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    | Convergency criteria checked for the last 2 iterations                       +
    | Convergence targets achieved. Performing final SCF iteration                 +
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     Total energy                               :      -152.41224082
     _______________________________________________________________
     Fermi energy                               :         0.33318841
     Kinetic energy                             :       152.06135667
     Coulomb energy                             :      -282.76363838
     Exchange energy                            :       -20.81913392
     Correlation energy                         :        -0.86315814
     DFT-D2 dispersion correction               :        -0.02766705
     
     DOS at Fermi energy (states/Ha/cell)       :        32.73559542
     
     Electron charges :
         core                                   :         8.00000000
         core leakage                           :         0.00261990
         valence                                :        16.00000000
         interstitial                           :         9.41909972
         charge in muffin-tin spheres :
                      atom     1     C          :         3.64065647
                      atom     2     C          :         3.64065647
                      atom     3     C          :         3.64979367
                      atom     4     C          :         3.64979367
         total charge in muffin-tins            :        14.58090028
         total charge                           :        24.00000000
     
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    + Self-consistent loop stopped                                                 +
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     STATE.OUT is written
     
    ********************************************************************************
    * Groundstate module stopped                                                   *
    ********************************************************************************
     
     Total time spent (seconds)                 :        47.78
    ================================================================================
    | EXCITING NEON stopped                                                        =
    ================================================================================
    """
    info_out_strain_2_str = info_out_strain_1_str.replace('-152.41224082', '-152.42761453')
    info_out_strain_3_str = info_out_strain_1_str.replace('-152.41224082', '-152.42597226')
    info_out_strain_4_str = info_out_strain_1_str.replace('-152.41224082', '-152.42367095')
    info_out_strain_oo_str = info_out_strain_1_str.replace('-152.41224082', '-152.41932764')

    nr_strain = 4
    for i_s in range(0, nr_strain):
        os.makedirs(os.path.dirname(tmp_path / f"rundir-{i_s + 1}/INFO.OUT"), exist_ok=True)

    os.makedirs(os.path.dirname(tmp_path / f"rundir-oo/INFO.OUT"), exist_ok=True)

    info_out_strain_1_file = tmp_path / "rundir-1/INFO.OUT"
    info_out_strain_1_file.write_text(info_out_strain_1_str)

    info_out_strain_2_file = tmp_path / "rundir-2/INFO.OUT"
    info_out_strain_2_file.write_text(info_out_strain_2_str)

    info_out_strain_3_file = tmp_path / "rundir-3/INFO.OUT"
    info_out_strain_3_file.write_text(info_out_strain_3_str)

    info_out_strain_4_file = tmp_path / "rundir-4/INFO.OUT"
    info_out_strain_4_file.write_text(info_out_strain_4_str)

    info_out_strain_oo_file = tmp_path / "rundir-oo/INFO.OUT"
    info_out_strain_oo_file.write_text(info_out_strain_oo_str)

    return MockFile(info_out_strain_1_file, info_out_strain_1_str), MockFile(info_out_strain_2_file, info_out_strain_2_str), \
           MockFile(info_out_strain_3_file, info_out_strain_3_str), MockFile(info_out_strain_4_file, info_out_strain_4_str), \
           MockFile(info_out_strain_oo_file, info_out_strain_oo_str)



def test_execute_elastic_strain(monkeypatch, info_out_mock, tmp_path):
    strain_values = [5, 6, 7, 8]
    def replace_run_exciting(rundir, *args, **kwargs):
        try:
            i = int(rundir.split('-')[-1])
            with open(info_out_mock[i-1].full_path, "w") as f:
                f.write(info_out_mock[i-1].string)

            with open(f"{tmp_path}/rundir-{i}/strain-{i}", "w") as f:
                f.write(str(strain_values[i-1]))
        except ValueError:
            with open(f"{tmp_path}/rundir-oo/strain-oo", "w") as f:
                f.write(str(20))

    monkeypatch.setattr(excitingscripts.execute.elastic_strain, 'run_exciting', replace_run_exciting)

    energy_strain_data_ref = np.array([[5, -152.4122408300],
                                        [6, -152.4276145300],
                                        [7, -152.4259722600],
                                        [8, -152.4236709500],
                                        [20, -152.4193276400]])

    assert np.allclose(execute_elastic_strain(root_directory=tmp_path), energy_strain_data_ref)

