#/usr/bin/env python
#PBS -l nodes=4:ppn=8
#PBS -l walltime=02:15:00

from ase import Atom, Atoms
from ase.io import read
from ase.constraints import FixAtoms
from ase.optimize.test import run_test
from gpaw import GPAW
from gpaw import Mixer
from gpaw.poisson import PoissonSolver

name = 'nanoparticle'

def get_atoms():
    atoms = Atoms([
        Atom('Pd', [5.078689759346383, 5.410678028467162, 4.000000000000000]),
        Atom('Pd', [7.522055777772603, 4.000000000000000, 4.000000000000000]),
        Atom('Pd', [7.522055777772603, 6.821356056934325, 4.000000000000000]),
        Atom('Pd', [6.707600438297196, 5.410678028467162, 6.303627574066606]),
        Atom('N',  [4.807604264052752, 5.728625577716107, 5.919407072553396]),
        Atom('H',  [4.000000000000000, 5.965167390141987, 6.490469524180266]),
    ])

    constraint = FixAtoms(mask=[a.symbol == 'Pd' for a in atoms])
    atoms.set_constraint(constraint)
    atoms.center(vacuum=4.0)
    atoms.set_pbc(False)
    return atoms

def get_calculator():
    calc = GPAW(gpts=(64, 64, 64), #h=0.18, gives 64x60x60
                mode='lcao', 
                basis='szp(dzp)',
                nbands=-5,
                xc='LDA',
                width=0.1,
                mixer=Mixer(beta=0.1, nmaxold=5, weight=50.0),
                poissonsolver=PoissonSolver(nn='M', relax='GS'),
                convergence={'energy': 1e-4, 'bands': -3},
                stencils=(3, 3),
                txt='nanoparticle.txt')
    return calc

run_test(get_atoms, get_calculator, name, fmax=0.05, steps=200)
