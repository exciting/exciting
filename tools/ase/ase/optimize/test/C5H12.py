#/usr/bin/env python
#PBS -l nodes=4:ppn=8
#PBS -l walltime=13:00:00
from ase import Atoms
from ase.optimize.test import run_test
from gpaw import GPAW
from gpaw import Mixer
from gpaw.poisson import PoissonSolver

name = 'C5H12'

def get_atoms():
    atoms = Atoms(symbols='C5H12',
                  pbc=[False, False, False],
                  cell=[
                      [ 16.83752497,   0.        ,   0.        ],
                      [  0.        ,  12.18645905,   0.        ],
                      [  0.        ,   0.        ,  11.83462179]
                  ],
                  positions=[
                      [  5.90380523,   5.65545388,   5.91569796],
                      [  7.15617518,   6.52907738,   5.91569796],
                      [  8.41815022,   5.66384716,   5.92196554],
                      [  9.68108996,   6.52891016,   5.91022362],
                      [ 10.93006206,   5.65545388,   5.91569796],
                      [  5.00000011,   6.30002353,   5.9163716 ],
                      [  5.88571848,   5.0122839 ,   6.82246859],
                      [  5.88625613,   5.01308931,   5.01214155],
                      [  7.14329342,   7.18115393,   6.81640316],
                      [  7.14551332,   7.17200869,   5.00879027],
                      [  8.41609966,   5.00661165,   5.02355167],
                      [  8.41971183,   5.0251482 ,   6.83462168],
                      [  9.69568096,   7.18645894,   6.8078633 ],
                      [  9.68914668,   7.16663649,   5.00000011],
                      [ 10.95518898,   5.02163182,   6.8289018 ],
                      [ 11.83752486,   6.29836826,   5.90274952],
                      [ 10.94464142,   5.00000011,   5.01802495]
                  ])
    return atoms

def get_calculator():
    calc = GPAW(h=0.2,
                mode = 'lcao',
                basis = 'szp(dzp)',
                mixer=Mixer(beta=0.1, nmaxold=5, weight=50.0),
                poissonsolver=PoissonSolver(nn='M', relax='GS'),
                txt='C5H12.txt')
    return calc

run_test(get_atoms, get_calculator, name + '-gpaw')
