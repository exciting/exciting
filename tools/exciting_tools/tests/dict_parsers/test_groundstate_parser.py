"""
Test ground state parsers
Execute tests from exciting_tools directory:
pytest --capture=tee-sys
"""

from excitingtools.exciting_dict_parsers.groundstate_parser import (
    parse_info_out,
    parse_linengy,
    parse_lo_recommendation,
)


def test_parse_info_out(tmp_path):
    """
    Note, this test will break if:
      * the parser keys change
      * the structure of the output changes
    """

    # Reference of the dictionary parsed with parse_info_out
    # Generated with print(json.dumps(info_out, sort_keys=True, indent=4))
    # Only retained first and last SCF keys
    info_ref = {
        "initialization": {
            "APW functions": 8,
            "Brillouin zone volume": 0.0734963595,
            "Effective Wigner radius, r_s": 3.55062021,
            "Exchange-correlation type": 100,
            "G-vector grid sizes": "36    36    36",
            "Lattice vectors (cartesian)": [15.0, 0.0, 0.0, 0.0, 15.0, 0.0, 0.0, 0.0, 15.0],
            "Maximum Hamiltonian size": 263,
            "Maximum number of plane-waves": 251,
            "Maximum |G+k| for APW functions": 1.66666667,
            "Maximum |G| for potential and density": 7.5,
            "Number of Bravais lattice symmetries": 48,
            "Number of crystal symmetries": 48,
            "Number of empty states": 5,
            "Polynomial order for pseudochg. density": 9,
            "R^MT_min * |G+k|_max (rgkmax)": 10.0,
            "Reciprocal lattice vectors (cartesian)": [
                0.4188790205,
                0.0,
                0.0,
                0.0,
                0.4188790205,
                0.0,
                0.0,
                0.0,
                0.4188790205,
            ],
            "Smearing scheme": "Gaussian",
            "Smearing width": 0.001,
            "Species 1": {
                "# of radial points in muffin-tin": 1000,
                "Atomic positions": {"Atom 1": "0.00000000  0.00000000  0.00000000"},
                "Species": "1 (Ar)",
                "Species symbol": "Ar",
                "atomic mass": 72820.7492,
                "electronic charge": 18.0,
                "muffin-tin radius": 6.0,
                "name": "argon",
                "nuclear charge": -18.0,
                "parameters loaded from": "Ar.xml",
            },
            "Species with R^MT_min": "1 (Ar)",
            "Spin treatment": "spin-unpolarised",
            "Total core charge": 10.0,
            "Total electronic charge": 18.0,
            "Total nuclear charge": -18.0,
            "Total number of G-vectors": 23871,
            "Total number of atoms per unit cell": 1,
            "Total number of k-points": 1,
            "Total number of local-orbitals": 12,
            "Total number of valence states": 10,
            "Total valence charge": 8.0,
            "Unit cell volume": 3375.0,
            "computing H and O matrix elements": 4,
            "inner part of muffin-tin": 2,
            "k-point grid": "1    1    1",
            "libxc; exchange": "Slater exchange; correlation",
            "mixing": "Using multisecant Broyden potential mixing",
            "potential and density": 4,
            "units": {"positions": "lattice"},
        },
        "scl": {
            "1": {
                "Core-electron kinetic energy": 0.0,
                "Correlation energy": -1.43085548,
                "Coulomb energy": -1029.02167746,
                "Coulomb potential energy": -796.81322609,
                "DOS at Fermi energy (states/Ha/cell)": 0.0,
                "Effective potential energy": -835.64023227,
                "Electron charges": "",
                "Electron-nuclear energy": -1208.12684923,
                "Estimated fundamental gap": 0.36071248,
                "Exchange energy": -27.93377198,
                "Fermi energy": -0.20111449,
                "Hartree energy": 205.65681157,
                "Kinetic energy": 530.56137212,
                "Madelung energy": -630.61506441,
                "Nuclear-nuclear energy": -26.5516398,
                "Sum of eigenvalues": -305.07886015,
                "Total energy": -527.82493279,
                "Wall time (seconds)": 1.05,
                "atom     1    Ar": 17.99816103,
                "charge in muffin-tin spheres": "",
                "core": 10.0,
                "core leakage": 0.0,
                "interstitial": 0.00183897,
                "total charge": 18.0,
                "total charge in muffin-tins": 17.99816103,
                "valence": 8.0,
                "xc potential energy": -38.82700618,
            },
            "12": {
                "Core-electron kinetic energy": 0.0,
                "Correlation energy": -1.4308435,
                "Coulomb energy": -1029.02642037,
                "Coulomb potential energy": -796.82023455,
                "DOS at Fermi energy (states/Ha/cell)": 0.0,
                "Effective potential energy": -835.64716936,
                "Electron charges": "",
                "Electron-nuclear energy": -1208.12932661,
                "Estimated fundamental gap": 0.36095838,
                "Exchange energy": -27.93372809,
                "Fermi energy": -0.20044598,
                "Hartree energy": 205.65454603,
                "Kinetic energy": 530.57303096,
                "Madelung energy": -630.6163031,
                "Nuclear-nuclear energy": -26.55163980,
                "Sum of eigenvalues": -305.0741384,
                "Total energy": -527.81796101,
                "atom     1    Ar": 17.99815963,
                "charge in muffin-tin spheres": "",
                "core": 10.0,
                "core leakage": 0.0,
                "interstitial": 0.00184037,
                "total charge": 18.0,
                "total charge in muffin-tins": 17.99815963,
                "valence": 8.0,
                "xc potential energy": -38.82693481,
            },
        },
    }

    file = tmp_path / "INFO.OUT"
    file.write_text(LDA_VWN_Ar_INFO_OUT)
    assert file.exists(), "INFO.OUT not written to tmp_path"

    info_out = parse_info_out(file.as_posix())

    assert info_out["initialization"] == info_ref["initialization"], "Initialization data consistent"
    assert len(info_out["scl"]) == 12, "expected 12 SCF steps"
    assert info_out["scl"]["1"] == info_ref["scl"]["1"], "SCF first iteration data consistent"
    assert info_out["scl"]["12"] == info_ref["scl"]["12"], "SCF last iteration data consistent"

    # Testing for the Structure Optimization Module
    # Only retained first and last SCF keys
    # Only retained first and last STR_OPT keys
    relax_info_ref = {
        "initialization": {
            "Lattice vectors (cartesian)": [10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 5.0],
            "Reciprocal lattice vectors (cartesian)": [
                0.6283185307,
                0.0,
                0.0,
                0.0,
                0.6283185307,
                0.0,
                0.0,
                0.0,
                1.2566370614,
            ],
            "Unit cell volume": 500.0,
            "Brillouin zone volume": 0.4961004269,
            "Species 1": {
                "Species symbol": "O",
                "Species": "1 (O)",
                "parameters loaded from": "O.xml",
                "name": "oxygen",
                "Atomic positions": {"Atom 1": "0.00000000  0.00000000  0.00000000"},
            },
            "Species 2": {
                "Species symbol": "H",
                "Species": "2 (H)",
                "parameters loaded from": "H.xml",
                "name": "hydrogen",
                "Atomic positions": {
                    "Atom 1": "2.50000000  0.50000000  0.00000000",
                    "Atom 2": "0.50000000  2.50000000  0.00000000",
                },
            },
            "Total number of atoms per unit cell": 3,
            "Spin treatment": "spin-unpolarised",
            "Number of Bravais lattice symmetries": 16,
            "Number of crystal symmetries": 4,
            "k-point grid": "1    1    1",
            "Total number of k-points": 1,
            "R^MT_min * |G+k|_max (rgkmax)": 3.0,
            "G-vector grid sizes": "54    54    27",
            "Total number of G-vectors": 34597,
            "Total nuclear charge": -10.0,
            "Total electronic charge": 10.0,
            "Number of empty states": 5,
            "Total number of valence states": 10,
            "Exchange-correlation type": 3,
            "Smearing scheme": "Gaussian",
            "Smearing width": 0.001,
            "units": {"positions": "cartesian"},
        },
        "scl": {
            "1": {
                "Total energy": -78.57899369,
                "Fermi energy": -0.17199772,
                "DOS at Fermi energy (states/Ha/cell)": 0.0,
                "Electron charge": ["Core", "leakage"],
                "Estimated fundamental gap": 0.10500197,
                "Wall time (seconds)": 2.06,
            },
            # ...
            "14": {
                "Total energy": -75.68495694,
                "Fermi energy": -0.09763242,
                "DOS at Fermi energy (states/Ha/cell)": 0.53544446,
                "Electron charge": ["Core", "leakage"],
            },
        },
        "str_opt": {
            0: {
                "Number of total scf iterations": 14,
                "Maximum force": [0.08541136, 0.0002],
                "Total energy": -75.68495694,
                "Atomic positions": {1: [0.0, 0.0, 0.0], 2: [2.5, 0.5, 0.0], 3: [0.5, 2.5, 0.0]},
                "Total atomic forces": {
                    1: [0.11374934, 0.11374934, -0.0],
                    2: [-0.07719269, -0.03655665, -0.0],
                    3: [-0.03655665, -0.07719269, 0.0],
                },
                "Center of mass": [0.16784752, 0.16784752, 0.0],
                "Total torque": [0.0, 0.0, -0.0],
            },
            # ...
            8: {
                "Number of total scf iterations": 13,
                "Maximum force": [5.82e-06, 0.0002],
                "Total energy": -75.78640119,
                "Atomic positions": {
                    1: [0.0, 0.0, 0.0],
                    2: [1.85302198, 9.80178632, 0.0],
                    3: [9.80178632, 1.85302198, 0.0],
                },
                "Total atomic forces": {
                    1: [7.04e-06, 7.04e-06, -0.0],
                    2: [-5.65e-06, -1.4e-06, 0.0],
                    3: [-1.4e-06, -5.65e-06, -0.0],
                },
                "Center of mass": [0.09258515, 0.09258515, 0.0],
                "Total torque": [-0.0, -0.0, -0.0],
            },
        },
    }

    file.write_text(LDA_PW_H2O_RELAX_INFO_OUT)
    assert file.exists(), "INFO.OUT not written to tmp_path"

    relax_info_out = parse_info_out(file.as_posix())

    assert len(relax_info_out["str_opt"]) == 9, "expected 9 relaxation steps"
    assert relax_info_out["str_opt"][0] == relax_info_ref["str_opt"][0], "Relaxation first step data consistent"
    assert relax_info_out["str_opt"][8] == relax_info_ref["str_opt"][8], "Relaxation last step data consistent"


LDA_VWN_Ar_INFO_OUT = """================================================================================
| EXCITING NITROGEN-14 started                                                 =
| version hash id: 1775bff4453c84689fb848894a9224f155377cfc                    =
|                                                                              =
|                                                                              =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
| Time (hh:mm:ss)   : 20:02:27                                                 =
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
     15.0000000000      0.0000000000      0.0000000000
      0.0000000000     15.0000000000      0.0000000000
      0.0000000000      0.0000000000     15.0000000000
 
 Reciprocal lattice vectors (cartesian) :
      0.4188790205      0.0000000000      0.0000000000
      0.0000000000      0.4188790205      0.0000000000
      0.0000000000      0.0000000000      0.4188790205
 
 Unit cell volume                           :    3375.0000000000
 Brillouin zone volume                      :       0.0734963595
 
 Species :    1 (Ar)
     parameters loaded from                 :    Ar.xml
     name                                   :    argon
     nuclear charge                         :     -18.00000000
     electronic charge                      :      18.00000000
     atomic mass                            :   72820.74920000
     muffin-tin radius                      :       6.00000000
     # of radial points in muffin-tin       :    1000
 
     atomic positions (lattice) :
       1 :   0.00000000  0.00000000  0.00000000
 
 Total number of atoms per unit cell        :       1
 
 Spin treatment                             :    spin-unpolarised
 
 Number of Bravais lattice symmetries       :      48
 Number of crystal symmetries               :      48
 
 k-point grid                               :       1    1    1
 Total number of k-points                   :       1
 k-point set is reduced with crystal symmetries
 
 R^MT_min * |G+k|_max (rgkmax)              :      10.00000000
 Species with R^MT_min                      :       1 (Ar)
 Maximum |G+k| for APW functions            :       1.66666667
 Maximum |G| for potential and density      :       7.50000000
 Polynomial order for pseudochg. density    :       9
 
 G-vector grid sizes                        :      36    36    36
 Total number of G-vectors                  :   23871
 
 Maximum angular momentum used for
     APW functions                          :       8
     computing H and O matrix elements      :       4
     potential and density                  :       4
     inner part of muffin-tin               :       2
 
 Total nuclear charge                       :     -18.00000000
 Total electronic charge                    :      18.00000000
 Total core charge                          :      10.00000000
 Total valence charge                       :       8.00000000
 
 Effective Wigner radius, r_s               :       3.55062021
 
 Number of empty states                     :       5
 Total number of valence states             :      10
 
 Maximum Hamiltonian size                   :     263
 Maximum number of plane-waves              :     251
 Total number of local-orbitals             :      12
 
 Exchange-correlation type                  :     100
     libxc; exchange: Slater exchange; correlation: Vosko, Wilk & Nusair (VWN5) (see libxc for references)
 
 Smearing scheme                            :    Gaussian
 Smearing width                             :       0.00100000
 
 Using multisecant Broyden potential mixing
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ Ending initialization                                                        +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
********************************************************************************
* Groundstate module started                                                   *
********************************************************************************
 Output level for this task is set to high
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ Self-consistent loop started                                                 +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Density and potential initialised from atomic data
 
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    1                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :      -527.82493279
 _______________________________________________________________
 Fermi energy                               :        -0.20111449
 Kinetic energy                             :       530.56137212
 Coulomb energy                             :     -1029.02167746
 Exchange energy                            :       -27.93377198
 Correlation energy                         :        -1.43085548
 Sum of eigenvalues                         :      -305.07886015
 Effective potential energy                 :      -835.64023227
 Coulomb potential energy                   :      -796.81322609
 xc potential energy                        :       -38.82700618
 Hartree energy                             :       205.65681157
 Electron-nuclear energy                    :     -1208.12684923
 Nuclear-nuclear energy                     :       -26.55163980
 Madelung energy                            :      -630.61506441
 Core-electron kinetic energy               :         0.00000000
 
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 
 Electron charges :
     core                                   :        10.00000000
     core leakage                           :         0.00000000
     valence                                :         8.00000000
     interstitial                           :         0.00183897
     charge in muffin-tin spheres :
                  atom     1    Ar          :        17.99816103
     total charge in muffin-tins            :        17.99816103
     total charge                           :        18.00000000
 
 Estimated fundamental gap                  :         0.36071248
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000
 
 Wall time (seconds)                        :         1.05
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    2                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :      -527.82194513
 _______________________________________________________________
 Fermi energy                               :        -0.20093816
 Kinetic energy                             :       530.56518416
 Coulomb energy                             :     -1029.02256672
 Exchange energy                            :       -27.93371266
 Correlation energy                         :        -1.43084992
 Sum of eigenvalues                         :      -305.07718688
 Effective potential energy                 :      -835.64237105
 Coulomb potential energy                   :      -796.81544990
 xc potential energy                        :       -38.82692115
 Hartree energy                             :       205.65547702
 Electron-nuclear energy                    :     -1208.12640394
 Nuclear-nuclear energy                     :       -26.55163980
 Madelung energy                            :      -630.61484177
 Core-electron kinetic energy               :         0.00000000
 
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 
 Electron charges :
     core                                   :        10.00000000
     core leakage                           :        -0.00000000
     valence                                :         8.00000000
     interstitial                           :         0.00183981
     charge in muffin-tin spheres :
                  atom     1    Ar          :        17.99816019
     total charge in muffin-tins            :        17.99816019
     total charge                           :        18.00000000
 
 Estimated fundamental gap                  :         0.36074458
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000
 
 Wall time (seconds)                        :         1.46
 
 RMS change in effective potential (target) :  0.587126E-02  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.298766E-02  ( 0.100000E-06)
 Charge distance                   (target) :  0.233904E-04  ( 0.100000E-04)
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    3                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :      -527.82022198
 _______________________________________________________________
 Fermi energy                               :        -0.20079854
 Kinetic energy                             :       530.56785510
 Coulomb energy                             :     -1029.02353844
 Exchange energy                            :       -27.93369195
 Correlation energy                         :        -1.43084669
 Sum of eigenvalues                         :      -305.07613298
 Effective potential energy                 :      -835.64398808
 Coulomb potential energy                   :      -796.81709801
 xc potential energy                        :       -38.82689007
 Hartree energy                             :       205.65480064
 Electron-nuclear energy                    :     -1208.12669929
 Nuclear-nuclear energy                     :       -26.55163980
 Madelung energy                            :      -630.61498944
 Core-electron kinetic energy               :         0.00000000
 
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 
 Electron charges :
     core                                   :        10.00000000
     core leakage                           :         0.00000000
     valence                                :         8.00000000
     interstitial                           :         0.00184026
     charge in muffin-tin spheres :
                  atom     1    Ar          :        17.99815974
     total charge in muffin-tins            :        17.99815974
     total charge                           :        18.00000000
 
 Estimated fundamental gap                  :         0.36078432
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000
 
 Wall time (seconds)                        :         1.82
 
 RMS change in effective potential (target) :  0.433810E-02  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.172315E-02  ( 0.100000E-06)
 Charge distance                   (target) :  0.150521E-04  ( 0.100000E-04)
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    4                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :      -527.82282966
 _______________________________________________________________
 Fermi energy                               :        -0.20048984
 Kinetic energy                             :       530.57007281
 Coulomb energy                             :     -1029.02813987
 Exchange energy                            :       -27.93391070
 Correlation energy                         :        -1.43085190
 Sum of eigenvalues                         :      -305.07638893
 Effective potential energy                 :      -835.64646174
 Coulomb potential energy                   :      -796.81927455
 xc potential energy                        :       -38.82718719
 Hartree energy                             :       205.65722553
 Electron-nuclear energy                    :     -1208.13372561
 Nuclear-nuclear energy                     :       -26.55163980
 Madelung energy                            :      -630.61850260
 Core-electron kinetic energy               :         0.00000000
 
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 
 Electron charges :
     core                                   :        10.00000000
     core leakage                           :        -0.00000000
     valence                                :         8.00000000
     interstitial                           :         0.00183901
     charge in muffin-tin spheres :
                  atom     1    Ar          :        17.99816099
     total charge in muffin-tins            :        17.99816099
     total charge                           :        18.00000000
 
 Estimated fundamental gap                  :         0.36103015
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000
 
 Wall time (seconds)                        :         2.15
 
 RMS change in effective potential (target) :  0.115090E-03  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.260768E-02  ( 0.100000E-06)
 Charge distance                   (target) :  0.157915E-04  ( 0.100000E-04)
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    5                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :      -527.81254947
 _______________________________________________________________
 Fermi energy                               :        -0.20040531
 Kinetic energy                             :       530.57593246
 Coulomb energy                             :     -1029.02412325
 Exchange energy                            :       -27.93352379
 Correlation energy                         :        -1.43083489
 Sum of eigenvalues                         :      -305.07152179
 Effective potential energy                 :      -835.64745425
 Coulomb potential energy                   :      -796.82080096
 xc potential energy                        :       -38.82665329
 Hartree energy                             :       205.65168249
 Electron-nuclear energy                    :     -1208.12416595
 Nuclear-nuclear energy                     :       -26.55163980
 Madelung energy                            :      -630.61372277
 Core-electron kinetic energy               :         0.00000000
 
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 
 Electron charges :
     core                                   :        10.00000000
     core leakage                           :         0.00000000
     valence                                :         8.00000000
     interstitial                           :         0.00184168
     charge in muffin-tin spheres :
                  atom     1    Ar          :        17.99815832
     total charge in muffin-tins            :        17.99815832
     total charge                           :        18.00000000
 
 Estimated fundamental gap                  :         0.36089040
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000
 
 Wall time (seconds)                        :         2.86
 
 RMS change in effective potential (target) :  0.144178E-03  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.102802E-01  ( 0.100000E-06)
 Charge distance                   (target) :  0.556814E-04  ( 0.100000E-04)
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    6                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :      -527.81742129
 _______________________________________________________________
 Fermi energy                               :        -0.20044051
 Kinetic energy                             :       530.57339029
 Coulomb energy                             :     -1029.02626113
 Exchange energy                            :       -27.93370794
 Correlation energy                         :        -1.43084251
 Sum of eigenvalues                         :      -305.07389872
 Effective potential energy                 :      -835.64728901
 Coulomb potential energy                   :      -796.82038212
 xc potential energy                        :       -38.82690689
 Hartree energy                             :       205.65423921
 Electron-nuclear energy                    :     -1208.12886054
 Nuclear-nuclear energy                     :       -26.55163980
 Madelung energy                            :      -630.61607007
 Core-electron kinetic energy               :         0.00000000
 
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 
 Electron charges :
     core                                   :        10.00000000
     core leakage                           :        -0.00000000
     valence                                :         8.00000000
     interstitial                           :         0.00184054
     charge in muffin-tin spheres :
                  atom     1    Ar          :        17.99815946
     total charge in muffin-tins            :        17.99815946
     total charge                           :        18.00000000
 
 Estimated fundamental gap                  :         0.36094948
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000
 
 Wall time (seconds)                        :         3.25
 
 RMS change in effective potential (target) :  0.117020E-04  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.487182E-02  ( 0.100000E-06)
 Charge distance                   (target) :  0.245233E-04  ( 0.100000E-04)
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    7                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :      -527.81797586
 _______________________________________________________________
 Fermi energy                               :        -0.20044613
 Kinetic energy                             :       530.57302094
 Coulomb energy                             :     -1029.02642463
 Exchange energy                            :       -27.93372865
 Correlation energy                         :        -1.43084353
 Sum of eigenvalues                         :      -305.07414492
 Effective potential energy                 :      -835.64716586
 Coulomb potential energy                   :      -796.82023028
 xc potential energy                        :       -38.82693558
 Hartree energy                             :       205.65455455
 Electron-nuclear energy                    :     -1208.12933939
 Nuclear-nuclear energy                     :       -26.55163980
 Madelung energy                            :      -630.61630949
 Core-electron kinetic energy               :         0.00000000
 
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 
 Electron charges :
     core                                   :        10.00000000
     core leakage                           :        -0.00000000
     valence                                :         8.00000000
     interstitial                           :         0.00184037
     charge in muffin-tin spheres :
                  atom     1    Ar          :        17.99815963
     total charge in muffin-tins            :        17.99815963
     total charge                           :        18.00000000
 
 Estimated fundamental gap                  :         0.36095863
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000
 
 Wall time (seconds)                        :         3.58
 
 RMS change in effective potential (target) :  0.323152E-06  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.554576E-03  ( 0.100000E-06)
 Charge distance                   (target) :  0.346569E-05  ( 0.100000E-04)
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    8                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :      -527.81796070
 _______________________________________________________________
 Fermi energy                               :        -0.20044598
 Kinetic energy                             :       530.57303118
 Coulomb energy                             :     -1029.02642030
 Exchange energy                            :       -27.93372808
 Correlation energy                         :        -1.43084350
 Sum of eigenvalues                         :      -305.07413828
 Effective potential energy                 :      -835.64716946
 Coulomb potential energy                   :      -796.82023466
 xc potential energy                        :       -38.82693479
 Hartree energy                             :       205.65454584
 Electron-nuclear energy                    :     -1208.12932634
 Nuclear-nuclear energy                     :       -26.55163980
 Madelung energy                            :      -630.61630297
 Core-electron kinetic energy               :         0.00000000
 
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 
 Electron charges :
     core                                   :        10.00000000
     core leakage                           :        -0.00000000
     valence                                :         8.00000000
     interstitial                           :         0.00184037
     charge in muffin-tin spheres :
                  atom     1    Ar          :        17.99815963
     total charge in muffin-tins            :        17.99815963
     total charge                           :        18.00000000
 
 Estimated fundamental gap                  :         0.36095838
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000
 
 Wall time (seconds)                        :         3.90
 
 RMS change in effective potential (target) :  0.740517E-08  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.151658E-04  ( 0.100000E-06)
 Charge distance                   (target) :  0.967001E-07  ( 0.100000E-04)
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    9                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :      -527.81796102
 _______________________________________________________________
 Fermi energy                               :        -0.20044598
 Kinetic energy                             :       530.57303095
 Coulomb energy                             :     -1029.02642037
 Exchange energy                            :       -27.93372809
 Correlation energy                         :        -1.43084350
 Sum of eigenvalues                         :      -305.07413841
 Effective potential energy                 :      -835.64716936
 Coulomb potential energy                   :      -796.82023455
 xc potential energy                        :       -38.82693481
 Hartree energy                             :       205.65454603
 Electron-nuclear energy                    :     -1208.12932661
 Nuclear-nuclear energy                     :       -26.55163980
 Madelung energy                            :      -630.61630310
 Core-electron kinetic energy               :         0.00000000
 
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 
 Electron charges :
     core                                   :        10.00000000
     core leakage                           :         0.00000000
     valence                                :         8.00000000
     interstitial                           :         0.00184037
     charge in muffin-tin spheres :
                  atom     1    Ar          :        17.99815963
     total charge in muffin-tins            :        17.99815963
     total charge                           :        18.00000000
 
 Estimated fundamental gap                  :         0.36095838
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000
 
 Wall time (seconds)                        :         4.23
 
 RMS change in effective potential (target) :  0.117161E-09  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.316730E-06  ( 0.100000E-06)
 Charge distance                   (target) :  0.220818E-08  ( 0.100000E-04)
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :   10                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :      -527.81796101
 _______________________________________________________________
 Fermi energy                               :        -0.20044598
 Kinetic energy                             :       530.57303096
 Coulomb energy                             :     -1029.02642037
 Exchange energy                            :       -27.93372809
 Correlation energy                         :        -1.43084350
 Sum of eigenvalues                         :      -305.07413840
 Effective potential energy                 :      -835.64716936
 Coulomb potential energy                   :      -796.82023455
 xc potential energy                        :       -38.82693481
 Hartree energy                             :       205.65454603
 Electron-nuclear energy                    :     -1208.12932661
 Nuclear-nuclear energy                     :       -26.55163980
 Madelung energy                            :      -630.61630310
 Core-electron kinetic energy               :         0.00000000
 
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 
 Electron charges :
     core                                   :        10.00000000
     core leakage                           :        -0.00000000
     valence                                :         8.00000000
     interstitial                           :         0.00184037
     charge in muffin-tin spheres :
                  atom     1    Ar          :        17.99815963
     total charge in muffin-tins            :        17.99815963
     total charge                           :        18.00000000
 
 Estimated fundamental gap                  :         0.36095838
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000
 
 Wall time (seconds)                        :         4.62
 
 RMS change in effective potential (target) :  0.250465E-10  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.419402E-08  ( 0.100000E-06)
 Charge distance                   (target) :  0.314056E-10  ( 0.100000E-04)
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :   11                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :      -527.81796101
 _______________________________________________________________
 Fermi energy                               :        -0.20044598
 Kinetic energy                             :       530.57303096
 Coulomb energy                             :     -1029.02642037
 Exchange energy                            :       -27.93372809
 Correlation energy                         :        -1.43084350
 Sum of eigenvalues                         :      -305.07413840
 Effective potential energy                 :      -835.64716936
 Coulomb potential energy                   :      -796.82023455
 xc potential energy                        :       -38.82693481
 Hartree energy                             :       205.65454603
 Electron-nuclear energy                    :     -1208.12932661
 Nuclear-nuclear energy                     :       -26.55163980
 Madelung energy                            :      -630.61630310
 Core-electron kinetic energy               :         0.00000000
 
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 
 Electron charges :
     core                                   :        10.00000000
     core leakage                           :        -0.00000000
     valence                                :         8.00000000
     interstitial                           :         0.00184037
     charge in muffin-tin spheres :
                  atom     1    Ar          :        17.99815963
     total charge in muffin-tins            :        17.99815963
     total charge                           :        18.00000000
 
 Estimated fundamental gap                  :         0.36095838
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000
 
 Wall time (seconds)                        :         4.95
 
 RMS change in effective potential (target) :  0.141030E-10  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.662567E-09  ( 0.100000E-06)
 Charge distance                   (target) :  0.430772E-11  ( 0.100000E-04)
                                                                                
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| Convergency criteria checked for the last 2 iterations                       +
| Convergence targets achieved. Performing final SCF iteration                 +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :      -527.81796101
 _______________________________________________________________
 Fermi energy                               :        -0.20044598
 Kinetic energy                             :       530.57303096
 Coulomb energy                             :     -1029.02642037
 Exchange energy                            :       -27.93372809
 Correlation energy                         :        -1.43084350
 Sum of eigenvalues                         :      -305.07413840
 Effective potential energy                 :      -835.64716936
 Coulomb potential energy                   :      -796.82023455
 xc potential energy                        :       -38.82693481
 Hartree energy                             :       205.65454603
 Electron-nuclear energy                    :     -1208.12932661
 Nuclear-nuclear energy                     :       -26.55163980
 Madelung energy                            :      -630.61630310
 Core-electron kinetic energy               :         0.00000000
 
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 
 Electron charges :
     core                                   :        10.00000000
     core leakage                           :         0.00000000
     valence                                :         8.00000000
     interstitial                           :         0.00184037
     charge in muffin-tin spheres :
                  atom     1    Ar          :        17.99815963
     total charge in muffin-tins            :        17.99815963
     total charge                           :        18.00000000
 
 Estimated fundamental gap                  :         0.36095838
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ Self-consistent loop stopped                                                 +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 STATE.OUT is written
 
********************************************************************************
* Groundstate module stopped                                                   *
********************************************************************************
 
--------------------------------------------------------------------------------
- Timings (seconds)                                                            -
--------------------------------------------------------------------------------
 Initialisation                             :         0.63
 Hamiltonian and overlap matrix set up      :         0.05
 First-variational secular equation         :         0.18
 Calculation of charge-density              :         0.14
 Calculation of potential                   :         2.97
 Muffin-tin manipulations                   :         1.02
 APW matching                               :         0.01
 Disk reads/writes                          :         0.37
 Mixing efforts                             :         0.02
 Solver of Dirac eqn.                       :         0.47
 Solver of rel. Schroedinger eqn.           :         0.48
 Total time spent in radial solvers         :         0.95
 
 Total time spent (seconds)                 :         4.18
================================================================================
| EXCITING NITROGEN-14 stopped                                                 =
================================================================================
"""

LDA_PW_H2O_RELAX_INFO_OUT = """================================================================================
| EXCITING NEON started                                                        =
| version hash id: 65d681e4e9a82d65818e11a67bb277686daac300                    =
|                                                                              =
| compiler: GNU Fortran (Debian 8.3.0-6) 8.3.0                                 =
|                                                                              =
|                                                                              =
| Date (DD-MM-YYYY) : 24-05-2024                                               =
| Time (hh:mm:ss)   : 11:43:08                                                 =
|                                                                              =
| All units are atomic (Hartree, Bohr, etc.)                                   =
================================================================================

********************************************************************************
* Structural optimisation starting from atomic densities                       *
********************************************************************************

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ Starting initialization                                                      +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 Lattice vectors (cartesian) :
     10.0000000000      0.0000000000      0.0000000000
      0.0000000000     10.0000000000      0.0000000000
      0.0000000000      0.0000000000      5.0000000000

 Reciprocal lattice vectors (cartesian) :
      0.6283185307      0.0000000000      0.0000000000
      0.0000000000      0.6283185307      0.0000000000
      0.0000000000      0.0000000000      1.2566370614

 Unit cell volume                           :     500.0000000000
 Brillouin zone volume                      :       0.4961004269

 Species :    1 (O)
     parameters loaded from                 :    O.xml
     name                                   :    oxygen

     atomic positions (cartesian) :
       1 :   0.00000000  0.00000000  0.00000000

 Species :    2 (H)
     parameters loaded from                 :    H.xml
     name                                   :    hydrogen

     atomic positions (cartesian) :
       1 :   2.50000000  0.50000000  0.00000000
       2 :   0.50000000  2.50000000  0.00000000

 Total number of atoms per unit cell        :       3

 Spin treatment                             :    spin-unpolarised

 Number of Bravais lattice symmetries       :      16
 Number of crystal symmetries               :       4

 k-point grid                               :       1    1    1
 Total number of k-points                   :       1
 k-point set is reduced with crystal symmetries

 R^MT_min * |G+k|_max (rgkmax)              :       3.00000000

 G-vector grid sizes                        :      54    54    27
 Total number of G-vectors                  :   34597

 Total nuclear charge                       :     -10.00000000
 Total electronic charge                    :      10.00000000

 Number of empty states                     :       5
 Total number of valence states             :      10

 Exchange-correlation type                  :       3
     Perdew-Wang, Phys. Rev. B 45, 13244 (1992)

 Smearing scheme                            :    Gaussian
 Smearing width                             :       0.00100000

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ Ending initialization                                                        +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

********************************************************************************
* Groundstate module started                                                   *
********************************************************************************
 Output level for this task is set to low

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ Self-consistent loop started                                                 +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Density and potential initialised from atomic data


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    1                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -78.57899369
 _______________________________________________________________
 Fermi energy                               :        -0.17199772
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 Electron charge: Core leakage              :         0.00015571

 Estimated fundamental gap                  :         0.10500197
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000

 Wall time (seconds)                        :         2.06

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    2                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -75.91027382
 _______________________________________________________________
 Fermi energy                               :        -0.12370214
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 Electron charge: Core leakage              :         0.00015663

 Estimated fundamental gap                  :         0.04196020
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000

 Wall time (seconds)                        :         3.36

 RMS change in effective potential (target) :  0.229615E-01  ( 0.100000E-05)
 Absolute change in total energy   (target) :   2.66872      ( 0.100000E-05)
 Charge distance                   (target) :  0.790693E-01  ( 0.100000E-04)
 Abs. change in max-nonIBS-force   (target) :  0.418141E-01  ( 0.500000E-04)

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    3                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -75.77423177
 _______________________________________________________________
 Fermi energy                               :        -0.11525374
 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 Electron charge: Core leakage              :         0.00015608

 Estimated fundamental gap                  :         0.03015303
        valence-band maximum at    1      0.0000  0.0000  0.0000
     conduction-band minimum at    1      0.0000  0.0000  0.0000

 Wall time (seconds)                        :         4.69

 RMS change in effective potential (target) :  0.165558E-01  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.136042      ( 0.100000E-05)
 Charge distance                   (target) :  0.468368E-02  ( 0.100000E-04)
 Abs. change in max-nonIBS-force   (target) :  0.274449E-01  ( 0.500000E-04)

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    4                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -75.71447271
 _______________________________________________________________
 Fermi energy                               :        -0.09763860
 DOS at Fermi energy (states/Ha/cell)       :         0.00576810
 Electron charge: Core leakage              :         0.00015463

 Wall time (seconds)                        :         7.26

 RMS change in effective potential (target) :  0.263203E-02  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.597591E-01  ( 0.100000E-05)
 Charge distance                   (target) :  0.815183E-02  ( 0.100000E-04)
 Abs. change in max-nonIBS-force   (target) :  0.130811E-01  ( 0.500000E-04)

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    5                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -75.64649264
 _______________________________________________________________
 Fermi energy                               :        -0.09727389
 DOS at Fermi energy (states/Ha/cell)       :         2.93046521
 Electron charge: Core leakage              :         0.00015446

 Wall time (seconds)                        :         8.65

 RMS change in effective potential (target) :  0.161339E-02  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.679801E-01  ( 0.100000E-05)
 Charge distance                   (target) :  0.455893E-02  ( 0.100000E-04)
 Abs. change in max-nonIBS-force   (target) :  0.604927E-02  ( 0.500000E-04)

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    6                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -75.68863676
 _______________________________________________________________
 Fermi energy                               :        -0.09762584
 DOS at Fermi energy (states/Ha/cell)       :         0.45697789
 Electron charge: Core leakage              :         0.00015447

 Wall time (seconds)                        :         9.77

 RMS change in effective potential (target) :  0.223804E-03  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.421441E-01  ( 0.100000E-05)
 Charge distance                   (target) :  0.179484E-02  ( 0.100000E-04)
 Abs. change in max-nonIBS-force   (target) :  0.439246E-02  ( 0.500000E-04)

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    7                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -75.68797134
 _______________________________________________________________
 Fermi energy                               :        -0.09767956
 DOS at Fermi energy (states/Ha/cell)       :         0.46423109
 Electron charge: Core leakage              :         0.00015448

 Wall time (seconds)                        :        10.93

 RMS change in effective potential (target) :  0.106267E-03  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.665420E-03  ( 0.100000E-05)
 Charge distance                   (target) :  0.159597E-03  ( 0.100000E-04)
 Abs. change in max-nonIBS-force   (target) :  0.125163E-02  ( 0.500000E-04)

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    8                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -75.68495355
 _______________________________________________________________
 Fermi energy                               :        -0.09763125
 DOS at Fermi energy (states/Ha/cell)       :         0.53598377
 Electron charge: Core leakage              :         0.00015448

 Wall time (seconds)                        :        12.06

 RMS change in effective potential (target) :  0.430409E-05  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.301779E-02  ( 0.100000E-05)
 Charge distance                   (target) :  0.946456E-04  ( 0.100000E-04)
 Abs. change in max-nonIBS-force   (target) :  0.773983E-04  ( 0.500000E-04)

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :    9                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -75.68494097
 _______________________________________________________________
 Fermi energy                               :        -0.09763220
 DOS at Fermi energy (states/Ha/cell)       :         0.53614420
 Electron charge: Core leakage              :         0.00015448

 Wall time (seconds)                        :        13.17

 RMS change in effective potential (target) :  0.733676E-06  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.125827E-04  ( 0.100000E-05)
 Charge distance                   (target) :  0.416156E-05  ( 0.100000E-04)
 Abs. change in max-nonIBS-force   (target) :  0.380398E-04  ( 0.500000E-04)

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :   10                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -75.68495103
 _______________________________________________________________
 Fermi energy                               :        -0.09763235
 DOS at Fermi energy (states/Ha/cell)       :         0.53551437
 Electron charge: Core leakage              :         0.00015448

 Wall time (seconds)                        :        14.48

 RMS change in effective potential (target) :  0.226552E-06  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.100633E-04  ( 0.100000E-05)
 Charge distance                   (target) :  0.555919E-06  ( 0.100000E-04)
 Abs. change in max-nonIBS-force   (target) :  0.319221E-05  ( 0.500000E-04)

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :   11                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -75.68495669
 _______________________________________________________________
 Fermi energy                               :        -0.09763241
 DOS at Fermi energy (states/Ha/cell)       :         0.53545372
 Electron charge: Core leakage              :         0.00015448

 Wall time (seconds)                        :        15.74

 RMS change in effective potential (target) :  0.111702E-07  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.565855E-05  ( 0.100000E-05)
 Charge distance                   (target) :  0.221398E-06  ( 0.100000E-04)
 Abs. change in max-nonIBS-force   (target) :  0.375697E-06  ( 0.500000E-04)

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :   12                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -75.68495693
 _______________________________________________________________
 Fermi energy                               :        -0.09763242
 DOS at Fermi energy (states/Ha/cell)       :         0.53544477
 Electron charge: Core leakage              :         0.00015448

 Wall time (seconds)                        :        17.00

 RMS change in effective potential (target) :  0.383891E-09  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.243353E-06  ( 0.100000E-05)
 Charge distance                   (target) :  0.917031E-08  ( 0.100000E-04)
 Abs. change in max-nonIBS-force   (target) :  0.107936E-07  ( 0.500000E-04)

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ SCF iteration number :   13                                                  +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -75.68495694
 _______________________________________________________________
 Fermi energy                               :        -0.09763242
 DOS at Fermi energy (states/Ha/cell)       :         0.53544434
 Electron charge: Core leakage              :         0.00015448

 Wall time (seconds)                        :        18.20

 RMS change in effective potential (target) :  0.119770E-09  ( 0.100000E-05)
 Absolute change in total energy   (target) :  0.483250E-08  ( 0.100000E-05)
 Charge distance                   (target) :  0.560836E-09  ( 0.100000E-04)
 Abs. change in max-nonIBS-force   (target) :  0.128180E-08  ( 0.500000E-04)
                                                                                
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| Convergency criteria checked for the last 2 iterations                       +
| Convergence targets achieved. Performing final SCF iteration                 +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -75.68495694
 _______________________________________________________________
 Fermi energy                               :        -0.09763242
 DOS at Fermi energy (states/Ha/cell)       :         0.53544446
 Electron charge: Core leakage              :         0.00015448

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ Self-consistent loop stopped                                                 +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 STATE.OUT is written

--------------------------------------------------------------------------------
- Writing atomic positions and forces                                          -
--------------------------------------------------------------------------------

 Atomic positions (cartesian) :
 atom     1   O  :     0.00000000    0.00000000    0.00000000
 atom     2   H  :     2.50000000    0.50000000    0.00000000
 atom     3   H  :     0.50000000    2.50000000    0.00000000

 Total atomic forces including IBS (cartesian) + constraints (cartesian) :
 atom     1   O  :     0.11374934    0.11374934   -0.00000000   T  T  T  
 atom     2   H  :    -0.07719269   -0.03655665   -0.00000000   F  F  T  
 atom     3   H  :    -0.03655665   -0.07719269    0.00000000   F  F  T  

 Atomic force components including IBS (cartesian) :
 atom     1   O  :    -0.07662130   -0.07662130   -0.00000000   HF force
                 :     0.15092973    0.15092973    0.00000000   core correction
                 :     0.03944092    0.03944092    0.00000000   IBS correction
 atom     2   H  :     0.01925792    0.05736338    0.00000000   HF force
                 :    -0.07546486   -0.07546486    0.00000000   core correction
                 :    -0.02098575   -0.01845517    0.00000000   IBS correction
 atom     3   H  :     0.05736338    0.01925792   -0.00000000   HF force
                 :    -0.07546486   -0.07546486    0.00000000   core correction
                 :    -0.01845517   -0.02098575    0.00000000   IBS correction

********************************************************************************
* Groundstate module stopped                                                   *
********************************************************************************

********************************************************************************
* Structure-optimization module started                                        *
********************************************************************************
 Output level for this task is set to normal

 Maximum displacement tau_BFGS   is         :         0.5000
 Maximum displacement tau_newton is         :         0.2000

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ Use L-BFGS-B method for optimizing atomic positions                          +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

--------------------------------------------------------------------------------
- Optimization step 0: Initialize optimization                                 -
--------------------------------------------------------------------------------
 Maximum force magnitude           (target) :         0.08541136  (  0.00020000)
 Total energy at this optimization step     :       -75.68495694

 Atomic positions at this step (cartesian) :
 atom     1   O  :     0.00000000    0.00000000    0.00000000
 atom     2   H  :     2.50000000    0.50000000    0.00000000
 atom     3   H  :     0.50000000    2.50000000    0.00000000

 Total atomic forces including IBS (cartesian) + constraints (cartesian) :
 atom     1   O  :     0.11374934    0.11374934   -0.00000000   T  T  T  
 atom     2   H  :    -0.07719269   -0.03655665   -0.00000000   F  F  T  
 atom     3   H  :    -0.03655665   -0.07719269    0.00000000   F  F  T  

 Center of mass  :     0.16784752    0.16784752    0.00000000  (cartesian)
 Total torque    :     0.00000000    0.00000000   -0.00000000  (cartesian)

--------------------------------------------------------------------------------
- Optimization step    1    (method = bfgs)                                    -
--------------------------------------------------------------------------------
 Number of investigated configurations      :         1
 Number of total scf iterations             :        19
 Maximum force magnitude           (target) :         0.08260579  (  0.00020000)
 Total energy at this optimization step     :       -75.69927763

 Atomic positions at this step (cartesian) :
 atom     1   O  :     0.00000000    0.00000000    0.00000000
 atom     2   H  :     2.42280731    0.46344335    0.00000000
 atom     3   H  :     0.46344335    2.42280731    0.00000000

 Total atomic forces including IBS (cartesian) + constraints (cartesian) :
 atom     1   O  :     0.11121513    0.11121513    0.00000000   T  T  T  
 atom     2   H  :    -0.07348652   -0.03772862    0.00000000   F  F  T  
 atom     3   H  :    -0.03772862   -0.07348652    0.00000000   F  F  T  

 Center of mass  :     0.16148333    0.16148333    0.00000000  (cartesian)
 Total torque    :     0.00000000    0.00000000    0.00000000  (cartesian)

 Time spent in this optimization step       :        21.12 seconds

--------------------------------------------------------------------------------
- Optimization step    2    (method = bfgs)                                    -
--------------------------------------------------------------------------------
 Number of investigated configurations      :         1
 Number of total scf iterations             :        16
 Maximum force magnitude           (target) :         0.04173782  (  0.00020000)
 Total energy at this optimization step     :       -75.77754003

 Atomic positions at this step (cartesian) :
 atom     1   O  :     0.00000000    0.00000000    0.00000000
 atom     2   H  :     2.00000000    0.00000000    0.00000000
 atom     3   H  :     0.00000000    2.00000000    0.00000000

 Total atomic forces including IBS (cartesian) + constraints (cartesian) :
 atom     1   O  :     0.05308190    0.05308190    0.00000000   T  T  T  
 atom     2   H  :    -0.03944834   -0.01363357    0.00000000   F  F  T  
 atom     3   H  :    -0.01363357   -0.03944834    0.00000000   F  F  T  

 Center of mass  :     0.11189834    0.11189834    0.00000000  (cartesian)
 Total torque    :     0.00000000    0.00000000    0.00000000  (cartesian)

 Time spent in this optimization step       :        19.26 seconds

--------------------------------------------------------------------------------
- Optimization step    3    (method = bfgs)                                    -
--------------------------------------------------------------------------------
 Number of investigated configurations      :         1
 Number of total scf iterations             :        18
 Maximum force magnitude           (target) :         0.03226320  (  0.00020000)
 Total energy at this optimization step     :       -75.78059284

 Atomic positions at this step (cartesian) :
 atom     1   O  :     0.00000000    0.00000000    0.00000000
 atom     2   H  :     1.96055166    9.98636643    0.00000000
 atom     3   H  :     9.98636643    1.96055166    0.00000000

 Total atomic forces including IBS (cartesian) + constraints (cartesian) :
 atom     1   O  :     0.04293570    0.04293570   -0.00000000   T  T  T  
 atom     2   H  :    -0.02918721   -0.01374849   -0.00000000   F  F  T  
 atom     3   H  :    -0.01374849   -0.02918721    0.00000000   F  F  T  

 Center of mass  :     0.10892846    0.10892846    0.00000000  (cartesian)
 Total torque    :     0.00000000    0.00000000    0.00000000  (cartesian)

 Time spent in this optimization step       :        19.79 seconds

--------------------------------------------------------------------------------
- Optimization step    4    (method = bfgs)                                    -
--------------------------------------------------------------------------------
 Number of investigated configurations      :         1
 Number of total scf iterations             :        15
 Maximum force magnitude           (target) :         0.01961100  (  0.00020000)
 Total energy at this optimization step     :       -75.78522038

 Atomic positions at this step (cartesian) :
 atom     1   O  :     0.00000000    0.00000000    0.00000000
 atom     2   H  :     1.82888799    9.88029496    0.00000000
 atom     3   H  :     9.88029496    1.82888799    0.00000000

 Total atomic forces including IBS (cartesian) + constraints (cartesian) :
 atom     1   O  :    -0.00717923   -0.00717923    0.00000000   T  T  T  
 atom     2   H  :     0.01698403   -0.00980480    0.00000000   F  F  T  
 atom     3   H  :    -0.00980480    0.01698403    0.00000000   F  F  T  

 Center of mass  :     0.09562737    0.09562737    0.00000000  (cartesian)
 Total torque    :     0.00000000    0.00000000    0.00000000  (cartesian)

 Time spent in this optimization step       :        15.92 seconds

--------------------------------------------------------------------------------
- Optimization step    5    (method = bfgs)                                    -
--------------------------------------------------------------------------------
 Number of investigated configurations      :         1
 Number of total scf iterations             :        14
 Maximum force magnitude           (target) :         0.00610696  (  0.00020000)
 Total energy at this optimization step     :       -75.78606196

 Atomic positions at this step (cartesian) :
 atom     1   O  :     0.00000000    0.00000000    0.00000000
 atom     2   H  :     1.85663579    9.85874660    0.00000000
 atom     3   H  :     9.85874660    1.85663579    0.00000000

 Total atomic forces including IBS (cartesian) + constraints (cartesian) :
 atom     1   O  :     0.00327447    0.00327447    0.00000000   T  T  T  
 atom     2   H  :     0.00235863   -0.00563310    0.00000000   F  F  T  
 atom     3   H  :    -0.00563310    0.00235863    0.00000000   F  F  T  

 Center of mass  :     0.09597422    0.09597422    0.00000000  (cartesian)
 Total torque    :     0.00000000    0.00000000    0.00000000  (cartesian)

 Time spent in this optimization step       :        14.27 seconds

--------------------------------------------------------------------------------
- Optimization step    6    (method = bfgs)                                    -
--------------------------------------------------------------------------------
 Number of investigated configurations      :         1
 Number of total scf iterations             :        14
 Maximum force magnitude           (target) :         0.00185118  (  0.00020000)
 Total energy at this optimization step     :       -75.78635670

 Atomic positions at this step (cartesian) :
 atom     1   O  :     0.00000000    0.00000000    0.00000000
 atom     2   H  :     1.85845115    9.82032564    0.00000000
 atom     3   H  :     9.82032564    1.85845115    0.00000000

 Total atomic forces including IBS (cartesian) + constraints (cartesian) :
 atom     1   O  :     0.00257432    0.00257432    0.00000000   T  T  T  
 atom     2   H  :    -0.00104916   -0.00152517    0.00000000   F  F  T  
 atom     3   H  :    -0.00152517   -0.00104916    0.00000000   F  F  T  

 Center of mass  :     0.09392617    0.09392617    0.00000000  (cartesian)
 Total torque    :     0.00000000    0.00000000    0.00000000  (cartesian)

 Time spent in this optimization step       :        14.14 seconds

--------------------------------------------------------------------------------
- Optimization step    7    (method = bfgs)                                    -
--------------------------------------------------------------------------------
 Number of investigated configurations      :         1
 Number of total scf iterations             :        14
 Maximum force magnitude           (target) :         0.00037962  (  0.00020000)
 Total energy at this optimization step     :       -75.78639967

 Atomic positions at this step (cartesian) :
 atom     1   O  :     0.00000000    0.00000000    0.00000000
 atom     2   H  :     1.85405751    9.80287979    0.00000000
 atom     3   H  :     9.80287979    1.85405751    0.00000000

 Total atomic forces including IBS (cartesian) + constraints (cartesian) :
 atom     1   O  :     0.00041566    0.00041566    0.00000000   T  T  T  
 atom     2   H  :    -0.00037772   -0.00003794   -0.00000000   F  F  T  
 atom     3   H  :    -0.00003794   -0.00037772    0.00000000   F  F  T  

 Center of mass  :     0.09270427    0.09270427    0.00000000  (cartesian)
 Total torque    :     0.00000000    0.00000000   -0.00000000  (cartesian)

 Time spent in this optimization step       :        14.76 seconds

--------------------------------------------------------------------------------
- Optimization step    8    (method = bfgs)                                    -
--------------------------------------------------------------------------------
 Number of investigated configurations      :         1
 Number of total scf iterations             :        13
 Maximum force magnitude           (target) :         0.00000582  (  0.00020000)
 Total energy at this optimization step     :       -75.78640119

 Atomic positions at this step (cartesian) :
 atom     1   O  :     0.00000000    0.00000000    0.00000000
 atom     2   H  :     1.85302198    9.80178632    0.00000000
 atom     3   H  :     9.80178632    1.85302198    0.00000000

 Total atomic forces including IBS (cartesian) + constraints (cartesian) :
 atom     1   O  :     0.00000704    0.00000704   -0.00000000   T  T  T  
 atom     2   H  :    -0.00000565   -0.00000140    0.00000000   F  F  T  
 atom     3   H  :    -0.00000140   -0.00000565   -0.00000000   F  F  T  

 Center of mass  :     0.09258515    0.09258515    0.00000000  (cartesian)
 Total torque    :    -0.00000000   -0.00000000   -0.00000000  (cartesian)

 Time spent in this optimization step       :        14.98 seconds

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ Force convergence target achieved                                            +
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Total energy                               :       -75.78640119
 _______________________________________________________________
 Fermi energy                               :        -0.07435889

 DOS at Fermi energy (states/Ha/cell)       :         0.00000000
 Electron charge: Core leakage              :         0.00015647

 Optimized atomic positions (cartesian) :
 atom     1   O  :     0.00000000    0.00000000    0.00000000
 atom     2   H  :     1.85302198    9.80178632    0.00000000
 atom     3   H  :     9.80178632    1.85302198    0.00000000

 Total atomic forces including IBS (cartesian) + constraints (cartesian) :
 atom     1   O  :     0.00000704    0.00000704   -0.00000000   T  T  T  
 atom     2   H  :    -0.00000565   -0.00000140    0.00000000   F  F  T  
 atom     3   H  :    -0.00000140   -0.00000565   -0.00000000   F  F  T  

 Atomic force components including IBS (cartesian) :
 atom     1   O  :    -0.24241492   -0.24241492    0.00000000   HF force
                 :     0.19280958    0.19280958    0.00000000   core correction
                 :     0.04961239    0.04961239   -0.00000000   IBS correction
 atom     2   H  :     0.12062339    0.12179153    0.00000000   HF force
                 :    -0.09640479   -0.09640479    0.00000000   core correction
                 :    -0.02422425   -0.02538814   -0.00000000   IBS correction
 atom     3   H  :     0.12179153    0.12062339    0.00000000   HF force
                 :    -0.09640479   -0.09640479    0.00000000   core correction
                 :    -0.02538814   -0.02422425    0.00000000   IBS correction

********************************************************************************
* Structure-optimization module stopped                                        *
********************************************************************************

 Total time spent (seconds)                 :       151.22

--------------------------------------------------------------------------------
- CAUTION! Warnings have been written in file WARNING.OUT !                    -
--------------------------------------------------------------------------------
================================================================================
| EXCITING NEON stopped                                                        =
================================================================================
"""


def test_parse_linengy(tmp_path):
    """
    Note, this test will break if:
      * the parser keys change
      * the structure of the output changes
    """

    # Reference of the dictionary parsed with parse_linengy
    # Generated with print(json.dumps(info_out, sort_keys=True, indent=4))

    linengy_ref = {
        "0": {
            "apw": [-2.1, -1.1, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15],
            "lo": [
                1.839014478,
                1.839014478,
                -0.7435250385,
                -0.7435250385,
                1.839014478,
                1.839014478,
                1.839014478,
                -1.775228393,
                -0.7435250385,
                -0.7435250385,
            ],
        },
        "1": {
            "apw": [-0.3, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15],
            "lo": [
                0.7475318132,
                0.7475318132,
                1.385527458,
                1.385527458,
                1.974765418,
                1.974765418,
                0.7475318132,
                0.7475318132,
                0.7475318132,
                -9.003002855,
                -9.003002855,
                -9.003002855,
                -9.003002855,
                -9.003002855,
                1.385527458,
                1.385527458,
                1.385527458,
                -6.768218591,
                -6.768218591,
                -6.768218591,
                -6.768218591,
                -6.768218591,
                1.974765418,
                1.974765418,
            ],
        },
    }

    file = tmp_path / "LINENGY.OUT"
    file.write_text(GGA_PBE_SOL_automatic_trial_energies_NaCl_LINENGY_OUT)
    assert file.exists(), "LINENGY.OUT not written to tmp_path"

    linengy = parse_linengy(file.as_posix())

    assert linengy["0"]["apw"] == linengy_ref["0"]["apw"], "First species apw data consistent"
    assert linengy["0"]["lo"] == linengy_ref["0"]["lo"], "First species lo data consistent"
    assert linengy["1"]["apw"] == linengy_ref["1"]["apw"], "Second species apw data consistent"
    assert linengy["1"]["lo"] == linengy_ref["1"]["lo"], "Second species lo data consistent"


GGA_PBE_SOL_automatic_trial_energies_NaCl_LINENGY_OUT = """Species :    1 (Na), atom :    1
 APW functions :
  l =  0, order =  1 :   -2.100000000    
  l =  1, order =  1 :   -1.100000000    
  l =  2, order =  1 :   0.1500000000    
  l =  2, order =  2 :   0.1500000000    
  l =  3, order =  1 :   0.1500000000    
  l =  3, order =  2 :   0.1500000000    
  l =  4, order =  1 :   0.1500000000    
  l =  4, order =  2 :   0.1500000000    
  l =  5, order =  1 :   0.1500000000    
  l =  5, order =  2 :   0.1500000000    
  l =  6, order =  1 :   0.1500000000    
  l =  6, order =  2 :   0.1500000000    
  l =  7, order =  1 :   0.1500000000    
  l =  7, order =  2 :   0.1500000000    
  l =  8, order =  1 :   0.1500000000    
  l =  8, order =  2 :   0.1500000000    
 local-orbital functions :
  l.o. =  1, l =  0, order =  1 :    1.839014478    
  l.o. =  1, l =  0, order =  2 :    1.839014478    
  l.o. =  2, l =  1, order =  1 :  -0.7435250385    
  l.o. =  2, l =  1, order =  2 :  -0.7435250385    
  l.o. =  3, l =  0, order =  1 :    1.839014478    
  l.o. =  3, l =  0, order =  2 :    1.839014478    
  l.o. =  4, l =  0, order =  1 :    1.839014478    
  l.o. =  4, l =  0, order =  2 :   -1.775228393    
  l.o. =  5, l =  1, order =  1 :  -0.7435250385    
  l.o. =  5, l =  1, order =  2 :  -0.7435250385    
 
Species :    2 (Cl), atom :    1
 APW functions :
  l =  0, order =  1 :  -0.3000000000    
  l =  1, order =  1 :   0.1500000000    
  l =  2, order =  1 :   0.1500000000    
  l =  3, order =  1 :   0.1500000000    
  l =  3, order =  2 :   0.1500000000    
  l =  4, order =  1 :   0.1500000000    
  l =  4, order =  2 :   0.1500000000    
  l =  5, order =  1 :   0.1500000000    
  l =  5, order =  2 :   0.1500000000    
  l =  6, order =  1 :   0.1500000000    
  l =  6, order =  2 :   0.1500000000    
  l =  7, order =  1 :   0.1500000000    
  l =  7, order =  2 :   0.1500000000    
  l =  8, order =  1 :   0.1500000000    
  l =  8, order =  2 :   0.1500000000    
 local-orbital functions :
  l.o. =  1, l =  0, order =  1 :   0.7475318132    
  l.o. =  1, l =  0, order =  2 :   0.7475318132    
  l.o. =  2, l =  1, order =  1 :    1.385527458    
  l.o. =  2, l =  1, order =  2 :    1.385527458    
  l.o. =  3, l =  2, order =  1 :    1.974765418    
  l.o. =  3, l =  2, order =  2 :    1.974765418    
  l.o. =  4, l =  0, order =  1 :   0.7475318132    
  l.o. =  4, l =  0, order =  2 :   0.7475318132    
  l.o. =  5, l =  0, order =  1 :   0.7475318132    
  l.o. =  5, l =  0, order =  2 :   -9.003002855    
  l.o. =  6, l =  0, order =  1 :   -9.003002855    
  l.o. =  6, l =  0, order =  2 :   -9.003002855    
  l.o. =  7, l =  0, order =  1 :   -9.003002855    
  l.o. =  7, l =  0, order =  2 :   -9.003002855    
  l.o. =  8, l =  1, order =  1 :    1.385527458    
  l.o. =  8, l =  1, order =  2 :    1.385527458    
  l.o. =  9, l =  1, order =  1 :    1.385527458    
  l.o. =  9, l =  1, order =  2 :   -6.768218591    
  l.o. = 10, l =  1, order =  1 :   -6.768218591    
  l.o. = 10, l =  1, order =  2 :   -6.768218591    
  l.o. = 11, l =  1, order =  1 :   -6.768218591    
  l.o. = 11, l =  1, order =  2 :   -6.768218591    
  l.o. = 12, l =  2, order =  1 :    1.974765418    
  l.o. = 12, l =  2, order =  2 :    1.974765418   """


def test_parse_lo_recommendation(tmp_path):
    """
    Note, this test will break if:
      * the parser keys change
      * the structure of the output changes
    """

    # Reference of the dictionary parsed with parse_lo_recommendation
    # Generated with print(json.dumps(info_out, sort_keys=True, indent=4))

    lo_recommendation_ref = {
        "n_species": 2,
        "n_l_channels": 4,
        "n_nodes": 21,
        "Na": {
            0: [
                [0.0, 1.0, -37.65990706181],
                [1.0, 2.0, -1.796282469826],
                [2.0, 3.0, 1.821425619362],
                [3.0, 4.0, 7.892575399809],
                [4.0, 5.0, 16.986165613637],
                [5.0, 6.0, 28.846080585981],
                [6.0, 7.0, 43.369227333017],
                [7.0, 8.0, 60.495491296112],
                [8.0, 9.0, 80.187724549122],
                [9.0, 10.0, 102.422176950914],
                [10.0, 11.0, 127.182990712735],
                [11.0, 12.0, 154.458866712915],
                [12.0, 13.0, 184.241158329361],
                [13.0, 14.0, 216.522955263927],
                [14.0, 15.0, 251.298756937161],
                [15.0, 16.0, 288.564346156074],
                [16.0, 17.0, 328.316602395923],
                [17.0, 18.0, 370.553203330272],
                [18.0, 19.0, 415.272297641907],
                [19.0, 20.0, 462.472263422762],
                [20.0, 21.0, 512.15160534914],
            ],
            1: [
                [0.0, 2.0, -0.764380938525],
                [1.0, 3.0, 2.150459202106],
                [2.0, 4.0, 7.68064315441],
                [3.0, 5.0, 15.944823175319],
                [4.0, 6.0, 26.826688781275],
                [5.0, 7.0, 40.283086869083],
                [6.0, 8.0, 56.289616891011],
                [7.0, 9.0, 74.830023777453],
                [8.0, 10.0, 95.891983053944],
                [9.0, 11.0, 119.465424264199],
                [10.0, 12.0, 145.542112894232],
                [11.0, 13.0, 174.115704751485],
                [12.0, 14.0, 205.181692970898],
                [13.0, 15.0, 238.737011824859],
                [14.0, 16.0, 274.779422707574],
                [15.0, 17.0, 313.306969990248],
                [16.0, 18.0, 354.317721540693],
                [17.0, 19.0, 397.809787603557],
                [18.0, 20.0, 443.781471675578],
                [19.0, 21.0, 492.231385973938],
                [20.0, 22.0, 543.158450851233],
            ],
            2: [
                [0.0, 3.0, 2.044659668292],
                [1.0, 4.0, 6.427853785834],
                [2.0, 5.0, 13.431600603035],
                [3.0, 6.0, 23.065707720667],
                [4.0, 7.0, 35.286723443086],
                [5.0, 8.0, 50.060639647089],
                [6.0, 9.0, 67.365856793742],
                [7.0, 10.0, 87.18903211822],
                [8.0, 11.0, 109.521849753452],
                [9.0, 12.0, 134.358761651081],
                [10.0, 13.0, 161.695405680071],
                [11.0, 14.0, 191.52770328581],
                [12.0, 15.0, 223.85164087452],
                [13.0, 16.0, 258.663517453378],
                [14.0, 17.0, 295.960290921283],
                [15.0, 18.0, 335.73972677672],
                [16.0, 19.0, 378.000263814603],
                [17.0, 20.0, 422.740724390955],
                [18.0, 21.0, 469.960046963641],
                [19.0, 22.0, 519.657164852307],
                [20.0, 23.0, 571.831020081621],
            ],
            3: [
                [0.0, 4.0, 3.879401688185],
                [1.0, 5.0, 9.964477889618],
                [2.0, 6.0, 18.47378344193],
                [3.0, 7.0, 29.505302500625],
                [4.0, 8.0, 43.071695950394],
                [5.0, 9.0, 59.165686791534],
                [6.0, 10.0, 77.776580159806],
                [7.0, 11.0, 98.89485358398],
                [8.0, 12.0, 122.513132271493],
                [9.0, 13.0, 148.626189112094],
                [10.0, 14.0, 177.230592656783],
                [11.0, 15.0, 208.324116670413],
                [12.0, 16.0, 241.905061644769],
                [13.0, 17.0, 277.971722458938],
                [14.0, 18.0, 316.522185899964],
                [15.0, 19.0, 357.554436071786],
                [16.0, 20.0, 401.066608041866],
                [17.0, 21.0, 447.057187991954],
                [18.0, 22.0, 495.525055433203],
                [19.0, 23.0, 546.469383913223],
                [20.0, 24.0, 599.889489594382],
            ],
        },
        "Cl": {
            0: [
                [0.0, 1.0, -100.974592263771],
                [1.0, 2.0, -8.943406336517],
                [2.0, 3.0, 0.803820215121],
                [3.0, 4.0, 11.265670674793],
                [4.0, 5.0, 27.995235804499],
                [5.0, 6.0, 50.139499480779],
                [6.0, 7.0, 77.425670325936],
                [7.0, 8.0, 109.722811813338],
                [8.0, 9.0, 146.943417860813],
                [9.0, 10.0, 189.030459361669],
                [10.0, 11.0, 235.945624307269],
                [11.0, 12.0, 287.660637659017],
                [12.0, 13.0, 344.155876726629],
                [13.0, 14.0, 405.415989381911],
                [14.0, 15.0, 471.42935525429],
                [15.0, 16.0, 542.186581729314],
                [16.0, 17.0, 617.680068606306],
                [17.0, 18.0, 697.90359882698],
                [18.0, 19.0, 782.852080659573],
                [19.0, 20.0, 872.521357718829],
                [20.0, 21.0, 966.908013775037],
            ],
            1: [
                [0.0, 2.0, -6.707821362791],
                [1.0, 3.0, 1.442002863349],
                [2.0, 4.0, 11.086329356142],
                [3.0, 5.0, 26.446761118127],
                [4.0, 6.0, 46.908960145565],
                [5.0, 7.0, 72.31449630133],
                [6.0, 8.0, 102.593748915541],
                [7.0, 9.0, 137.708872214526],
                [8.0, 10.0, 177.628868217528],
                [9.0, 11.0, 222.335815545055],
                [10.0, 12.0, 271.812774967963],
                [11.0, 13.0, 326.047911895729],
                [12.0, 14.0, 385.030908226461],
                [13.0, 15.0, 448.753672165588],
                [14.0, 16.0, 517.209705744261],
                [15.0, 17.0, 590.393819776316],
                [16.0, 18.0, 668.301879780751],
                [17.0, 19.0, 750.930453348811],
                [18.0, 20.0, 838.276595190814],
                [19.0, 21.0, 930.337672332757],
                [20.0, 22.0, 1027.111295649087],
            ],
            2: [
                [0.0, 3.0, 2.030626399034],
                [1.0, 4.0, 9.540398069814],
                [2.0, 5.0, 22.521806870793],
                [3.0, 6.0, 40.580505794574],
                [4.0, 7.0, 63.608967480673],
                [5.0, 8.0, 91.521851721479],
                [6.0, 9.0, 124.28110055295],
                [7.0, 10.0, 161.847878698612],
                [8.0, 11.0, 204.202486065259],
                [9.0, 12.0, 251.327905737439],
                [10.0, 13.0, 303.212514505293],
                [11.0, 14.0, 359.847427715492],
                [12.0, 15.0, 421.225227525549],
                [13.0, 16.0, 487.339835071361],
                [14.0, 17.0, 558.185801503932],
                [15.0, 18.0, 633.758349847869],
                [16.0, 19.0, 714.053292909356],
                [17.0, 20.0, 799.067051302979],
                [18.0, 21.0, 888.796639877965],
                [19.0, 22.0, 983.239587157865],
                [20.0, 23.0, 1082.393826570133],
            ],
            3: [
                [0.0, 4.0, 5.877668530234],
                [1.0, 5.0, 16.734445897029],
                [2.0, 6.0, 32.521343858776],
                [3.0, 7.0, 53.248636736501],
                [4.0, 8.0, 78.839760672149],
                [5.0, 9.0, 109.265751512702],
                [6.0, 10.0, 144.502223464379],
                [7.0, 11.0, 184.526132806013],
                [8.0, 12.0, 229.322756357998],
                [9.0, 13.0, 278.877755682806],
                [10.0, 14.0, 333.181582371119],
                [11.0, 15.0, 392.226089160902],
                [12.0, 16.0, 456.00524371215],
                [13.0, 17.0, 524.514261786364],
                [14.0, 18.0, 597.749258861014],
                [15.0, 19.0, 675.706956131225],
                [16.0, 20.0, 758.384407267418],
                [17.0, 21.0, 845.778892453188],
                [18.0, 22.0, 937.887883257667],
                [19.0, 23.0, 1034.709065214155],
                [20.0, 24.0, 1136.240377687629],
            ],
        },
    }

    file = tmp_path / "LO_RECOMMENDATION.OUT"
    file.write_text(GGA_PBE_SOL_automatic_trial_energies_NaCl_LO_RECOMMENDATION_OUT)
    assert file.exists(), "LO_RECOMMENDATION.OUT not written to tmp_path"

    lo_recommendation = parse_lo_recommendation(file.as_posix())

    assert lo_recommendation["Na"] == lo_recommendation_ref["Na"], "First species data consistent"
    assert lo_recommendation["Cl"] == lo_recommendation_ref["Cl"], "Second species data consistent"


GGA_PBE_SOL_automatic_trial_energies_NaCl_LO_RECOMMENDATION_OUT = """# Recommended linearization energies computet with Wigner-Seitz rules.
 --------------------------------------------------------------------
 #  n_species:  2
 # n_l-channels:  4
 # n_nodes: 21
 
 # species: Na, l :  0
 # nodes   n        trial energy
       0   1    -37.659907061810
       1   2     -1.796282469826
       2   3      1.821425619362
       3   4      7.892575399809
       4   5     16.986165613637
       5   6     28.846080585981
       6   7     43.369227333017
       7   8     60.495491296112
       8   9     80.187724549122
       9  10    102.422176950914
      10  11    127.182990712735
      11  12    154.458866712915
      12  13    184.241158329361
      13  14    216.522955263927
      14  15    251.298756937161
      15  16    288.564346156074
      16  17    328.316602395923
      17  18    370.553203330272
      18  19    415.272297641907
      19  20    462.472263422762
      20  21    512.151605349140
 
 # species: Na, l :  1
 # nodes   n        trial energy
       0   2     -0.764380938525
       1   3      2.150459202106
       2   4      7.680643154410
       3   5     15.944823175319
       4   6     26.826688781275
       5   7     40.283086869083
       6   8     56.289616891011
       7   9     74.830023777453
       8  10     95.891983053944
       9  11    119.465424264199
      10  12    145.542112894232
      11  13    174.115704751485
      12  14    205.181692970898
      13  15    238.737011824859
      14  16    274.779422707574
      15  17    313.306969990248
      16  18    354.317721540693
      17  19    397.809787603557
      18  20    443.781471675578
      19  21    492.231385973938
      20  22    543.158450851233
 
 # species: Na, l :  2
 # nodes   n        trial energy
       0   3      2.044659668292
       1   4      6.427853785834
       2   5     13.431600603035
       3   6     23.065707720667
       4   7     35.286723443086
       5   8     50.060639647089
       6   9     67.365856793742
       7  10     87.189032118220
       8  11    109.521849753452
       9  12    134.358761651081
      10  13    161.695405680071
      11  14    191.527703285810
      12  15    223.851640874520
      13  16    258.663517453378
      14  17    295.960290921283
      15  18    335.739726776720
      16  19    378.000263814603
      17  20    422.740724390955
      18  21    469.960046963641
      19  22    519.657164852307
      20  23    571.831020081621
 
 # species: Na, l :  3
 # nodes   n        trial energy
       0   4      3.879401688185
       1   5      9.964477889618
       2   6     18.473783441930
       3   7     29.505302500625
       4   8     43.071695950394
       5   9     59.165686791534
       6  10     77.776580159806
       7  11     98.894853583980
       8  12    122.513132271493
       9  13    148.626189112094
      10  14    177.230592656783
      11  15    208.324116670413
      12  16    241.905061644769
      13  17    277.971722458938
      14  18    316.522185899964
      15  19    357.554436071786
      16  20    401.066608041866
      17  21    447.057187991954
      18  22    495.525055433203
      19  23    546.469383913223
      20  24    599.889489594382
 
 # species: Cl, l :  0
 # nodes   n        trial energy
       0   1   -100.974592263771
       1   2     -8.943406336517
       2   3      0.803820215121
       3   4     11.265670674793
       4   5     27.995235804499
       5   6     50.139499480779
       6   7     77.425670325936
       7   8    109.722811813338
       8   9    146.943417860813
       9  10    189.030459361669
      10  11    235.945624307269
      11  12    287.660637659017
      12  13    344.155876726629
      13  14    405.415989381911
      14  15    471.429355254290
      15  16    542.186581729314
      16  17    617.680068606306
      17  18    697.903598826980
      18  19    782.852080659573
      19  20    872.521357718829
      20  21    966.908013775037
 
 # species: Cl, l :  1
 # nodes   n        trial energy
       0   2     -6.707821362791
       1   3      1.442002863349
       2   4     11.086329356142
       3   5     26.446761118127
       4   6     46.908960145565
       5   7     72.314496301330
       6   8    102.593748915541
       7   9    137.708872214526
       8  10    177.628868217528
       9  11    222.335815545055
      10  12    271.812774967963
      11  13    326.047911895729
      12  14    385.030908226461
      13  15    448.753672165588
      14  16    517.209705744261
      15  17    590.393819776316
      16  18    668.301879780751
      17  19    750.930453348811
      18  20    838.276595190814
      19  21    930.337672332757
      20  22   1027.111295649087
 
 # species: Cl, l :  2
 # nodes   n        trial energy
       0   3      2.030626399034
       1   4      9.540398069814
       2   5     22.521806870793
       3   6     40.580505794574
       4   7     63.608967480673
       5   8     91.521851721479
       6   9    124.281100552950
       7  10    161.847878698612
       8  11    204.202486065259
       9  12    251.327905737439
      10  13    303.212514505293
      11  14    359.847427715492
      12  15    421.225227525549
      13  16    487.339835071361
      14  17    558.185801503932
      15  18    633.758349847869
      16  19    714.053292909356
      17  20    799.067051302979
      18  21    888.796639877965
      19  22    983.239587157865
      20  23   1082.393826570133
 
 # species: Cl, l :  3
 # nodes   n        trial energy
       0   4      5.877668530234
       1   5     16.734445897029
       2   6     32.521343858776
       3   7     53.248636736501
       4   8     78.839760672149
       5   9    109.265751512702
       6  10    144.502223464379
       7  11    184.526132806013
       8  12    229.322756357998
       9  13    278.877755682806
      10  14    333.181582371119
      11  15    392.226089160902
      12  16    456.005243712150
      13  17    524.514261786364
      14  18    597.749258861014
      15  19    675.706956131225
      16  20    758.384407267418
      17  21    845.778892453188
      18  22    937.887883257667
      19  23   1034.709065214155
      20  24   1136.240377687629
"""
