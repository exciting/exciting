#!/usr/bin/python

"""
Run VASP tests to ensure that relaxation with the VASP calculator works.
This is conditional on the existence of the VASP_COMMAND or VASP_SCRIPT
environment variables.

"""

from ase.test import NotAvailable
import os

vcmd = os.getenv('VASP_COMMAND')
vscr = os.getenv('VASP_SCRIPT')
if vcmd == None and vscr == None:
    raise NotAvailable('Neither VASP_COMMAND nor VASP_SCRIPT defined')

import numpy as np
from ase import io
# QuasiNewton nowadays is an alias for BFGSLineSearch, which is
# broken. Use BFGS instead.
from ase.optimize import BFGS as QuasiNewton
from ase.lattice import bulk
from ase.calculators.vasp import Vasp

# -- Perform Volume relaxation within Vasp
def vasp_vol_relax():
    Al = bulk('Al', 'fcc', a=4.5, cubic=True)
    calc = Vasp(xc='LDA', isif=7, nsw=5,
                ibrion=1, ediffg=-1e-3, lwave=False, lcharg=False)
    calc.calculate(Al)

    # Explicitly parse atomic position output file from Vasp
    CONTCAR_Al = io.read('CONTCAR', format='vasp')

    print 'Stress after relaxation:\n', calc.read_stress()

    print 'Al cell post relaxation from calc:\n', calc.get_atoms().get_cell()
    print 'Al cell post relaxation from atoms:\n', Al.get_cell()
    print 'Al cell post relaxation from CONTCAR:\n', CONTCAR_Al.get_cell()

    # All the cells should be the same.
    assert (calc.get_atoms().get_cell() == CONTCAR_Al.get_cell()).all()
    assert (Al.get_cell() == CONTCAR_Al.get_cell()).all()

    return Al

# -- Perform Volume relaxation using ASE with Vasp as force/stress calculator
def ase_vol_relax():
    Al = bulk('Al', 'fcc', a=4.5, cubic=True)
    calc = Vasp(xc='LDA')
    Al.set_calculator(calc)

    from ase.constraints import StrainFilter
    sf = StrainFilter(Al)
    qn = QuasiNewton(sf, logfile='relaxation.log')
    qn.run(fmax=0.1, steps=5)

    print 'Stress:\n', calc.read_stress()
    print 'Al post ASE volume relaxation\n', calc.get_atoms().get_cell()

    return Al

# Test function for comparing two cells
def cells_almost_equal(cellA, cellB, tol=0.01):
    return  (np.abs(cellA - cellB) < tol).all()

# Correct LDA relaxed cell
a_rel = 4.18
LDA_cell = np.diag([a_rel, a_rel, a_rel])

Al_vasp = vasp_vol_relax()
Al_ase = ase_vol_relax()

assert cells_almost_equal(LDA_cell, Al_vasp.get_cell())
assert cells_almost_equal(LDA_cell, Al_ase.get_cell())

# Cleanup
Al_ase.get_calculator().clean()
