"""Effective medium theory potential."""

from math import sqrt, exp, log, pi

import numpy as np
import sys

from ase.data import chemical_symbols
from ase.units import Bohr
from ase.calculators.neighborlist import NeighborList


parameters = {
    #      E0     s0    V0     eta2    kappa   lambda  n0
    #      eV     bohr  eV     bohr^-1 bohr^-1 bohr^-1 bohr^-3
    'Al': (-3.28, 3.00, 1.493, 1.240,  2.000,  1.169,  0.00700),
    'Cu': (-3.51, 2.67, 2.476, 1.652,  2.740,  1.906,  0.00910),
    'Ag': (-2.96, 3.01, 2.132, 1.652,  2.790,  1.892,  0.00547),
    'Au': (-3.80, 3.00, 2.321, 1.674,  2.873,  2.182,  0.00703),
    'Ni': (-4.44, 2.60, 3.673, 1.669,  2.757,  1.948,  0.01030),
    'Pd': (-3.90, 2.87, 2.773, 1.818,  3.107,  2.155,  0.00688),
    'Pt': (-5.85, 2.90, 4.067, 1.812,  3.145,  2.192,  0.00802),
    # extra parameters - just for fun ...
    'H':  (-3.21, 1.31, 0.132, 2.652,  2.790,  3.892,  0.00547),
    'C':  (-3.50, 1.81, 0.332, 1.652,  2.790,  1.892,  0.01322),
    'N':  (-5.10, 1.88, 0.132, 1.652,  2.790,  1.892,  0.01222),
    'O':  (-4.60, 1.95, 0.332, 1.652,  2.790,  1.892,  0.00850)}

beta = 1.809     # (16 * pi / 3)**(1.0 / 3) / 2**0.5,
                 # but preserve historical rounding


class EMT:
    disabled = False  # Set to True to disable (asap does this).
    
    def __init__(self, fakestress=False):
        self.energy = None
        self.name = 'EMT'
        self.version = '1.0'
        # fakestress is needed to fake some stress value for the testsuite
        # in order to test the filter functionality.
        self.fakestress = fakestress
        if self.disabled:
            print >> sys.stderr, """
            ase.EMT has been disabled by Asap.  Most likely, you
            intended to use Asap's EMT calculator, but accidentally
            imported ase's EMT calculator after Asap's.  This could
            happen if your script contains the lines

              from asap3 import *
              from ase.calculators.emt import EMT
            Swap the two lines to solve the problem.

            (or 'from ase import *' in older versions of ASE.) Swap 
            the two lines to solve the problem.

            In the UNLIKELY event that you actually wanted to use
            ase.calculators.emt.EMT although asap3 is loaded into memory,
            please reactivate it with the command
              ase.calculators.emt.EMT.disabled = False
            """
            raise RuntimeError('ase.EMT has been disabled.  ' +
                               'See message printed above.')
        
    def get_name(self):
        return self.name

    def get_version(self):
        return self.version

    def get_spin_polarized(self):
        return False
    
    def initialize(self, atoms):
        self.par = {}
        self.rc = 0.0
        self.numbers = atoms.get_atomic_numbers()
        maxseq = max(par[1] for par in parameters.values()) * Bohr
        rc = self.rc = beta * maxseq * 0.5 * (sqrt(3) + sqrt(4))
        rr = rc * 2 * sqrt(4) / (sqrt(3) + sqrt(4))
        self.acut = np.log(9999.0) / (rr - rc)
        for Z in self.numbers:
            if Z not in self.par:
                p = parameters[chemical_symbols[Z]]
                s0 = p[1] * Bohr
                eta2 = p[3] / Bohr
                kappa = p[4] / Bohr
                x = eta2 * beta * s0
                gamma1 = 0.0
                gamma2 = 0.0
                for i, n in enumerate([12, 6, 24]):
                    r = s0 * beta * sqrt(i + 1)
                    x = n / (12 * (1.0 + exp(self.acut * (r - rc))))
                    gamma1 += x * exp(-eta2 * (r - beta * s0))
                    gamma2 += x * exp(-kappa / beta * (r - beta * s0))

                self.par[Z] = {'E0': p[0],
                               's0': s0,
                               'V0': p[2],
                               'eta2': eta2,
                               'kappa': kappa,
                               'lambda': p[5] / Bohr,
                               'n0': p[6] / Bohr**3,
                               'rc': rc,
                               'gamma1': gamma1,
                               'gamma2': gamma2}
                #if rc + 0.5 > self.rc:
                #    self.rc = rc + 0.5

        self.ksi = {}
        for s1, p1 in self.par.items():
            self.ksi[s1] = {}
            for s2, p2 in self.par.items():
                #self.ksi[s1][s2] = (p2['n0'] / p1['n0'] *
                #                    exp(eta1 * (p1['s0'] - p2['s0'])))
                self.ksi[s1][s2] = p2['n0'] / p1['n0']
                
        self.forces = np.empty((len(atoms), 3))
        self.sigma1 = np.empty(len(atoms))
        self.deds = np.empty(len(atoms))
                    
        self.nl = NeighborList([0.5 * self.rc + 0.25] * len(atoms),
                               self_interaction=False)

    def update(self, atoms):
        if (self.energy is None or
            len(self.numbers) != len(atoms) or
            (self.numbers != atoms.get_atomic_numbers()).any()):
            self.initialize(atoms)
            self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any() or
              (self.pbc != atoms.get_pbc()).any() or
              (self.cell != atoms.get_cell()).any()):
            self.calculate(atoms)

    def calculation_required(self, atoms, quantities):
        if len(quantities) == 0:
            return False

        return (self.energy is None or
                len(self.numbers) != len(atoms) or
                (self.numbers != atoms.get_atomic_numbers()).any() or
                (self.positions != atoms.get_positions()).any() or
                (self.pbc != atoms.get_pbc()).any() or
                (self.cell != atoms.get_cell()).any())
                
    def get_number_of_iterations(self):
        return 0

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy

    def get_numeric_forces(self, atoms):
        self.update(atoms)
        p = atoms.positions
        p0 = p.copy()
        forces = np.empty_like(p)
        eps = 0.0001
        for a in range(len(p)):
            for c in range(3):
                p[a, c] += eps
                self.calculate(atoms)
                de = self.energy
                p[a, c] -= 2 * eps
                self.calculate(atoms)
                de -= self.energy
                p[a, c] += eps
                forces[a, c] = -de / (2 * eps)
        p[:] = p0
        return forces

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()
    
    def get_stress(self, atoms):
        if self.fakestress:
            return np.zeros((6))
        else:
            raise NotImplementedError
    
    def calculate(self, atoms):
        self.positions = atoms.get_positions().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()
        
        self.nl.update(atoms)
        
        self.energy = 0.0
        self.sigma1[:] = 0.0
        self.forces[:] = 0.0
        
        natoms = len(atoms)

        for a1 in range(natoms):
            Z1 = self.numbers[a1]
            p1 = self.par[Z1]
            ksi = self.ksi[Z1]
            neighbors, offsets = self.nl.get_neighbors(a1)
            offsets = np.dot(offsets, atoms.cell)
            for a2, offset in zip(neighbors, offsets):
                d = self.positions[a2] + offset - self.positions[a1]
                r = sqrt(np.dot(d, d))
                if r < self.rc + 0.5:
                    Z2 = self.numbers[a2]
                    p2 = self.par[Z2]
                    self.interact1(a1, a2, d, r, p1, p2, ksi[Z2])
                                
        for a in range(natoms):
            Z = self.numbers[a]
            p = self.par[Z]
            try:
                ds = -log(self.sigma1[a] / 12) / (beta * p['eta2'])
            except (OverflowError, ValueError):
                self.deds[a] = 0.0
                self.energy -= p['E0']
                continue
            x = p['lambda'] * ds
            y = exp(-x)
            z = 6 * p['V0'] * exp(-p['kappa'] * ds)
            self.deds[a] = ((x * y * p['E0'] * p['lambda'] + p['kappa'] * z) /
                            (self.sigma1[a] * beta * p['eta2']))
            e = p['E0'] * ((1 + x) * y - 1) + z
            self.energy += p['E0'] * ((1 + x) * y - 1) + z

        for a1 in range(natoms):
            Z1 = self.numbers[a1]
            p1 = self.par[Z1]
            ksi = self.ksi[Z1]
            neighbors, offsets = self.nl.get_neighbors(a1)
            offsets = np.dot(offsets, atoms.cell)
            for a2, offset in zip(neighbors, offsets):
                d = self.positions[a2] + offset - self.positions[a1]
                r = sqrt(np.dot(d, d))
                if r < self.rc + 0.5:
                    Z2 = self.numbers[a2]
                    p2 = self.par[Z2]
                    self.interact2(a1, a2, d, r, p1, p2, ksi[Z2])

    def interact1(self, a1, a2, d, r, p1, p2, ksi):
        x = exp(self.acut * (r - self.rc))
        theta = 1.0 / (1.0 + x)
        y1 = (0.5 * p1['V0'] * exp(-p2['kappa'] * (r / beta - p2['s0'])) *
              ksi / p1['gamma2'] * theta)
        y2 = (0.5 * p2['V0'] * exp(-p1['kappa'] * (r / beta - p1['s0'])) /
              ksi / p2['gamma2'] * theta)
        self.energy -= y1 + y2
        f = ((y1 * p2['kappa'] + y2 * p1['kappa']) / beta +
             (y1 + y2) * self.acut * theta * x) * d / r
        self.forces[a1] += f
        self.forces[a2] -= f
        self.sigma1[a1] += (exp(-p2['eta2'] * (r - beta * p2['s0'])) *
                            ksi * theta / p1['gamma1'])
        self.sigma1[a2] += (exp(-p1['eta2'] * (r - beta * p1['s0'])) /
                            ksi * theta / p2['gamma1'])

    def interact2(self, a1, a2, d, r, p1, p2, ksi):
        x = exp(self.acut * (r - self.rc))
        theta = 1.0 / (1.0 + x)
        y1 = (exp(-p2['eta2'] * (r - beta * p2['s0'])) *
              ksi / p1['gamma1'] * theta * self.deds[a1])
        y2 = (exp(-p1['eta2'] * (r - beta * p1['s0'])) /
              ksi / p2['gamma1'] * theta * self.deds[a2])
        f = ((y1 * p2['eta2'] + y2 * p1['eta2']) +
             (y1 + y2) * self.acut * theta * x) * d / r
        self.forces[a1] -= f
        self.forces[a2] += f

    def set_atoms(self,*args,**kwargs):
        'empty function for compatibility with other calculators and tests'
        pass

