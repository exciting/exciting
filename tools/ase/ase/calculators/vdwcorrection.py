"""van der Waals correction schemes for DFT"""

import numpy as np
from ase.units import Bohr, Hartree
from gpaw.mpi import rank

# dipole polarizabilities and C6 values from 
# X. Chu and A. Dalgarno, J. Chem. Phys. 129 (2004) 4083
# atomic units, a_0^3
vdWDB_Chu04jcp = {
    # Element: [alpha, C6]; units [Bohr^3, Hartree * Bohr^6]
    'H'  : [4.5, 6.5], # [exact, Tkatchenko PRL]
    'He' : [1.38, 1.42],
    'Li' : [164, 1392],
    'Be' : [38, 227],
    'B'  : [21, 99.5],
    'C'  : [12, 46.6],
    'N'  : [7.4, 24.2],
    'O'  : [5.4, 15.6],
    'F'  : [3.8, 9.52],
    'Ne' : [2.67, 6.20],
    'Na' : [163, 1518],
    'Mg' : [71, 626],
    'Al' : [60, 528],
    'Si' : [37, 305],
    'Cl' : [15, 94.6],
    'Ar' : [11.1, 64.2],
    'Ca' : [160, 2163],
    'Fe' : [56, 482],
    'Br' : [20, 162],
    'Kr' : [16.7, 130],
    'Sr' : [199, 3175],
    'I'  : [35, 385],
}

# C6 values and vdW radii from 
# S. Grimme, J Comput Chem 27 (2006) 1787-1799
vdWDB_Grimme06jcc = {
    # Element: [C6, R0]; units [J nm^6 mol^{-1}, Angstrom]
    'H'  : [0.14, 1.001],
    'He' : [0.08, 1.012],
    'Li' : [1.61, 0.825],
    'Be' : [1.61, 1.408],
    'B'  : [3.13, 1.485],
    'C'  : [1.75, 1.452],
    'N'  : [1.23, 1.397],
    'O'  : [0.70, 1.342],
    'F'  : [0.75, 1.287],
    'Ne' : [0.63, 1.243],
    'Na' : [5.71, 1.144],
    'Mg' : [5.71, 1.364],
    'Al' : [10.79, 1.639],
    'Si' : [9.23, 1.716],
    'P'  : [7.84, 1.705],
    'S'  : [5.57, 1.683],
    'Cl' : [5.07, 1.639],
    'Ar' : [4.61, 1.595],
    'K'  : [10.80, 1.485],
    'Ca' : [10.80, 1.474],
    'Sc' : [10.80, 1.562],
    'Ti' : [10.80, 1.562],
    'V'  : [10.80, 1.562],
    'Cr'  : [10.80, 1.562],
    'Mn'  : [10.80, 1.562],
    'Fe'  : [10.80, 1.562],
    'Co'  : [10.80, 1.562],
    'Ni'  : [10.80, 1.562],
    'Cu'  : [10.80, 1.562],
    'Zn' : [10.80, 1.562],
    'Ga' : [16.99, 1.650],
    'Ge' : [17.10, 1.727],
    'As' : [16.37, 1.760],
    'Se' : [12.64, 1.771],
    'Br' : [12.47, 1.749],
    'Kr' : [12.01, 1.727],
    'Rb' : [24.67, 1.628],
    'Sr' : [24.67, 1.606],
    'Y-Cd' : [24.67, 1.639],
    'In' : [37.32, 1.672],
    'Sn' : [38.71, 1.804],
    'Sb' : [38.44, 1.881],
    'Te' : [31.74, 1.892],
    'I'  : [31.50, 1.892],
    'Xe' : [29.99, 1.881],
    }

class vdWTkatchenko09prl:
    """vdW correction after Tkatchenko and Scheffler PRL 102 (2009) 073005.

    hirshfeld: the Hirshfeld partitioning object
    calculator: the calculator to get the PBE energy
    """
    def __init__(self,                  
                 hirshfeld=None, vdwradii=None, calculator=None,
                 Rmax = 10, # maximal radius for periodic calculations
                 vdWDB_alphaC6 = vdWDB_Chu04jcp, # 
                 ):
        self.hirshfeld = hirshfeld
        if calculator is None:
            self.calculator = self.hirshfeld.get_calculator()
        else:
            self.calculator = calculator
        self.vdwradii = vdwradii
        self.vdWDB_alphaC6 = vdWDB_alphaC6
        self.Rmax = Rmax
        self.atoms = None

        self.sR = 0.94
        self.d = 20

    def update(self, atoms=None):
        if atoms is None:
            atoms = self.calculator.get_atoms()
        if (self.atoms and 
            (self.atoms.get_positions() == atoms.get_positions()).all()):
            return
        self.energy = self.calculator.get_potential_energy(atoms)
        self.forces = self.calculator.get_forces(atoms)
        self.atoms = atoms.copy()

        if self.vdwradii is not None:
            # external vdW radii
            vdwradii = self.vdwradii
            assert(len(atoms) == len(vdwradii))
        else:
            vdwradii = []
            for atom in atoms:
                self.vdwradii.append(vdWDB_Grimme06jcc[atom.symbol][1])
 
        if self.hirshfeld == None:
            volume_ratios = [1.] * len(atoms)
        elif hasattr(self.hirshfeld,'__len__'): # a list
            assert(len(atoms) == len(self.hirshfeld))
            volume_ratios = self.hirshfeld
        else: # sould be an object
            self.hirshfeld.initialize()
            volume_ratios = self.hirshfeld.get_effective_volume_ratios()

        # correction for effective C6
        na = len(atoms)
        C6eff_a = np.empty((na))
        alpha_a = np.empty((na))
        R0eff_a = np.empty((na))
        for a, atom in enumerate(atoms):
            # free atom values
            alpha_a[a], C6eff_a[a] = self.vdWDB_alphaC6[atom.symbol]
            # correction for effective C6
            C6eff_a[a] *= Hartree * volume_ratios[a]**2 * Bohr**6
            R0eff_a[a] = vdwradii[a] * volume_ratios[a]**(1./3.)
        C6eff_aa = np.empty((na, na))
        for a in range(na):
            for b in range(a, na):
                C6eff_aa[a, b] = (2 * C6eff_a[a] * C6eff_a[b] /
                                  (alpha_a[b] / alpha_a[a] * C6eff_a[a] +
                                   alpha_a[a] / alpha_a[b] * C6eff_a[b]   ))
                C6eff_aa[b, a] = C6eff_aa[a, b]

        # PBC
        pbc_c = atoms.get_pbc()
        cell_c = atoms.get_cell()
        Rcell_c = np.sqrt(np.sum(cell_c**2, axis=1))
        ncells_c = np.ceil(np.where(pbc_c, 1. + self.Rmax / Rcell_c, 1))
        ncells_c = np.array(ncells_c, dtype=int)

        positions = atoms.get_positions()
        EvdW = 0.0
        forces = 0. * self.forces
        # loop over all atoms in the cell
        for ia, posa in enumerate(positions):
            # loop over all atoms in the cell (and neighbour cells for PBC)
            for ib, posb in enumerate(positions):
                # loops over neighbour cells
                for ix in range(-ncells_c[0] + 1, ncells_c[0]):
                    for iy in range(-ncells_c[1] + 1, ncells_c[1]):
                        for iz in range(-ncells_c[2] + 1, ncells_c[2]):
                            i_c = np.array([ix, iy, iz])
                            diff = posb + np.dot(i_c, cell_c) - posa
                            r2 = np.dot(diff, diff)
                            r6 = r2**3
                            r = np.sqrt(r2)
                            if r > 1.e-10 and r < self.Rmax:
                                Edamp, Fdamp = self.damping(r, 
                                                            R0eff_a[ia],
                                                            R0eff_a[ib],
                                                            d=self.d, 
                                                            sR=self.sR)
                                EvdW -= (Edamp *
                                         C6eff_aa[ia, ib] / r6 )
                                # we neglect the C6eff contribution to the
                                # forces
                                forces[ia] -= ((Fdamp - 6 * Edamp / r) *
                                               C6eff_aa[ia, ib] / r6 *
                                               diff / r                 )
        self.energy += EvdW / 2. # double counting
        self.forces += forces / 2. # double counting
        
    def damping(self, RAB, R0A, R0B,
                d = 20,   # steepness of the step function
                sR = 0.94 # for PBE
                ):
        """Damping factor.

        Standard values for d and sR as given in 
        Tkatchenko and Scheffler PRL 102 (2009) 073005."""
        scale = 1.0 / (sR * (R0A + R0B))
        x = RAB * scale
        chi = np.exp(-d * (x - 1.0))
        return 1.0 / (1.0 + chi), d * scale * chi / (1.0 + chi)**2
 
    def get_potential_energy(self, atoms=None):
        self.update(atoms)
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces

    def get_stress(self, atoms):
        return np.zeros((3, 3))
