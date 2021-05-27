"""
Test ground state parsers
Execute tests from exciting_tools directory:
pytest --capture=tee-sys
"""

from excitingtools.parser.groundStateParser import parse_info_out


def test_parse_info_out():
    """
    Note, this test will break if:
      * the reference data location changes
      * the parser keys change
      * the structure of the output changes
    """

    # Reference of the dictionary parsed with parse_info_out
    # Generated with print(json.dumps(info_out, sort_keys=True, indent=4))
    # Only retained first and last SCF keys
    info_ref = {
        "initialization": {
            "# of radial points in muffin-tin": "1000",
            "1": "0.00000000  0.00000000  0.00000000",
            "APW functions": "8",
            "Brillouin zone volume": "0.0734963595",
            "Effective Wigner radius, r_s": "3.55062021",
            "Exchange-correlation type": "100",
            "G-vector grid sizes": "36    36    36",
            "Lattice vectors (cartesian)": [
                "15.0000000000",
                "0.0000000000",
                "0.0000000000",
                "0.0000000000",
                "15.0000000000",
                "0.0000000000",
                "0.0000000000",
                "0.0000000000",
                "15.0000000000"
            ],
            "Maximum Hamiltonian size": "263",
            "Maximum number of plane-waves": "251",
            "Maximum |G| for potential and density": "7.50000000",
            "Number of Bravais lattice symmetries": "48",
            "Number of crystal symmetries": "48",
            "Number of empty states": "5",
            "Polynomial order for pseudochg. density": "9",
            "Reciprocal lattice vectors (cartesian)": [
                "0.4188790205",
                "0.0000000000",
                "0.0000000000",
                "0.0000000000",
                "0.4188790205",
                "0.0000000000",
                "0.0000000000",
                "0.0000000000",
                "0.4188790205"
            ],
            "Smearing scheme": "Gaussian",
            "Smearing width": "0.00100000",
            "Species": "1 (Ar)",
            "Species with R^MT_min": "1 (Ar)",
            "Spin treatment": "spin-unpolarised",
            "Total core charge": "10.00000000",
            "Total electronic charge": "18.00000000",
            "Total nuclear charge": "-18.00000000",
            "Total number of G-vectors": "23871",
            "Total number of atoms per unit cell": "1",
            "Total number of k-points": "1",
            "Total number of local-orbitals": "12",
            "Total number of valence states": "10",
            "Total valence charge": "8.00000000",
            "Unit cell volume": "3375.0000000000",
            "atomic mass": "72820.74920000",
            "atomic positions (lattice)": "",
            "computing H and O matrix elements": "4",
            "electronic charge": "18.00000000",
            "inner part of muffin-tin": "2",
            "k-point grid": "1    1    1",
            "libxc; exchange": "Slater exchange; correlation",
            "mixing": "Using multisecant Broyden potential mixing",
            "muffin-tin radius": "6.00000000",
            "name": "argon",
            "nuclear charge": "-18.00000000",
            "parameters loaded from": "Ar.xml",
            "potential and density": "4"
        },
        "scl": {
            "1": {
                "Core-electron kinetic energy": "0.00000000",
                "Correlation energy": "-1.43085548",
                "Coulomb energy": "-1029.02167746",
                "Coulomb potential energy": "-796.81322609",
                "DOS at Fermi energy (states/Ha/cell)": "0.00000000",
                "Effective potential energy": "-835.64023227",
                "Electron charges": "",
                "Electron-nuclear energy": "-1208.12684923",
                "Estimated fundamental gap": "0.36071248",
                "Exchange energy": "-27.93377198",
                "Fermi energy": "-0.20111449",
                "Hartree energy": "205.65681157",
                "Kinetic energy": "530.56137212",
                "Madelung energy": "-630.61506441",
                "Nuclear-nuclear energy": "-26.55163980",
                "Sum of eigenvalues": "-305.07886015",
                "Total energy": "-527.82493279",
                "Wall time (seconds)": "1.05",
                "atom     1    Ar": "17.99816103",
                "charge in muffin-tin spheres": "",
                "core": "10.00000000",
                "core leakage": "0.00000000",
                "interstitial": "0.00183897",
                "total charge": "18.00000000",
                "total charge in muffin-tins": "17.99816103",
                "valence": "8.00000000",
                "xc potential energy": "-38.82700618"
            },
            "11": {
                "Core-electron kinetic energy": "0.00000000",
                "Correlation energy": "-1.43084350",
                "Coulomb energy": "-1029.02642037",
                "Coulomb potential energy": "-796.82023455",
                "DOS at Fermi energy (states/Ha/cell)": "0.00000000",
                "Effective potential energy": "-835.64716936",
                "Electron charges": "",
                "Electron-nuclear energy": "-1208.12932661",
                "Estimated fundamental gap": "0.36095838",
                "Exchange energy": "-27.93372809",
                "Fermi energy": "-0.20044598",
                "Hartree energy": "205.65454603",
                "Kinetic energy": "530.57303096",
                "Madelung energy": "-630.61630310",
                "Nuclear-nuclear energy": "-26.55163980",
                "Sum of eigenvalues": "-305.07413840",
                "Total energy": "-527.81796101",
                "Wall time (seconds)": "4.95",
                "atom     1    Ar": "17.99815963",
                "charge in muffin-tin spheres": "",
                "core": "10.00000000",
                "core leakage": "0.00000000",
                "interstitial": "0.00184037",
                "total charge": "18.00000000",
                "total charge in muffin-tins": "17.99815963",
                "valence": "8.00000000",
                "xc potential energy": "-38.82693481"
            }
        }
    }

    info_out = parse_info_out("../../test/test_farm/groundstate-GGA_PBE-Ar/ref/INFO.OUT")

    assert info_out['initialization'] == info_ref['initialization'], "Initialization data consistent"
    assert info_out['scl']['1'] == info_ref['scl']['1'], "SCF first iteration data consistent"
    assert info_out['scl']['11'] == info_ref['scl']['11'], "SCF last iteration data consistent"
