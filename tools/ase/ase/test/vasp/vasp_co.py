#!/usr/bin/python

"""
Run some VASP tests to ensure that the VASP calculator works. This
is conditional on the existence of the VASP_COMMAND or VASP_SCRIPT
environment variables

"""

from ase.test import NotAvailable
import os

vcmd = os.getenv('VASP_COMMAND')
vscr = os.getenv('VASP_SCRIPT')
if vcmd == None and vscr == None:
    raise NotAvailable('Neither VASP_COMMAND nor VASP_SCRIPT defined')

from ase import Atoms
from ase.calculators.vasp import Vasp
import numpy as np

def array_almost_equal(a1, a2, tol=np.finfo(type(1.0)).eps):
    """Replacement for old numpy.testing.utils.array_almost_equal."""
    return (np.abs(a1 - a2) < tol).all()

d = 1.14
co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)],
              pbc=True)
co.center(vacuum=5.)

calc = Vasp(
            xc = 'PBE',
            prec = 'Low',
            algo = 'Fast',
            ismear= 0,
            sigma = 1.,
            istart = 0,
            lwave = False,
            lcharg = False)

co.set_calculator(calc)
en = co.get_potential_energy()
assert abs(en + 14.918933) < 1e-4

# Secondly, check that restart from the previously created VASP output works

calc2 = Vasp(restart=True)
co2 = calc2.get_atoms()

# Need tolerance of 1e-14 because VASP itself changes coordinates
# slightly between reading POSCAR and writing CONTCAR even if no ionic
# steps are made.
assert array_almost_equal(co.positions, co2.positions, 1e-14)

assert en - co2.get_potential_energy() == 0.
assert array_almost_equal(calc.get_stress(co), calc2.get_stress(co2))
assert array_almost_equal(calc.get_forces(co), calc2.get_forces(co2))
assert array_almost_equal(calc.get_eigenvalues(), calc2.get_eigenvalues())
assert calc.get_number_of_bands() == calc2.get_number_of_bands()
assert calc.get_xc_functional() == calc2.get_xc_functional()

# Cleanup
calc.clean()
