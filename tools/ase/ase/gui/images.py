from math import sqrt

import numpy as np

from ase.data import covalent_radii
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write, string2index
from ase.constraints import FixAtoms
from ase.gui.defaults import read_defaults
from ase.quaternions import Quaternion

class Images:
    def __init__(self, images=None):

        if images is not None:
            self.initialize(images)
    
    def initialize(self, images, filenames=None, init_magmom=False):
        
        self.natoms = len(images[0])
        self.nimages = len(images)
        if hasattr(images[0], 'get_shapes'):
            self.shapes = images[0].get_shapes()
            self.Q = []
        else:
            self.shapes = None

        if filenames is None:
            filenames = [None] * self.nimages
        self.filenames = filenames
        self.P = np.empty((self.nimages, self.natoms, 3))
        self.V = np.empty((self.nimages, self.natoms, 3))
        self.E = np.empty(self.nimages)
        self.K = np.empty(self.nimages)
        self.F = np.empty((self.nimages, self.natoms, 3))
        self.M = np.empty((self.nimages, self.natoms))
        self.T = np.empty((self.nimages, self.natoms), int)
        self.A = np.empty((self.nimages, 3, 3))
        self.Z = images[0].get_atomic_numbers()
        self.pbc = images[0].get_pbc()
        self.covalent_radii = covalent_radii
        config = read_defaults()
        if config['covalent_radii'] is not None:
            for data in config['covalent_radii']:
                self.covalent_radii[data[0]] = data[1]
        warning = False
        for i, atoms in enumerate(images):
            natomsi = len(atoms)
            if (natomsi != self.natoms or
                (atoms.get_atomic_numbers() != self.Z).any()):
                raise RuntimeError('Can not handle different images with ' +
                                   'different numbers of atoms or different ' +
                                   'kinds of atoms!')
            self.P[i] = atoms.get_positions()
            self.V[i] = atoms.get_velocities()

            if hasattr(self, 'Q'):
                for q in atoms.get_quaternions():
                     self.Q.append(Quaternion(q))

            self.A[i] = atoms.get_cell()
            if (atoms.get_pbc() != self.pbc).any():
                warning = True
            try:
                self.E[i] = atoms.get_potential_energy()
            except RuntimeError:
                self.E[i] = np.nan
            self.K[i] = atoms.get_kinetic_energy()
            try:
                self.F[i] = atoms.get_forces(apply_constraint=False)
            except RuntimeError:
                self.F[i] = np.nan
            try:
                if init_magmom:
                    self.M[i] = atoms.get_initial_magnetic_moments()
                else:
                    self.M[i] = atoms.get_magnetic_moments()
            except (RuntimeError, AttributeError):
                self.M[i] = atoms.get_initial_magnetic_moments()
                
            # added support for tags
            try:
                self.T[i] = atoms.get_tags()
            except RuntimeError:
                self.T[i] = 0
                

        if warning:
            print('WARNING: Not all images have the same bondary conditions!')
            
        self.selected = np.zeros(self.natoms, bool)
        self.selected_ordered  = []
        self.atoms_to_rotate_0 = np.zeros(self.natoms, bool)
        self.visible = np.ones(self.natoms, bool)
        self.nselected = 0
        self.set_dynamic(constraints = images[0].constraints)
        self.repeat = np.ones(3, int)
        self.set_radii(config['radii_scale'])
        
    def prepare_new_atoms(self):
        "Marks that the next call to append_atoms should clear the images."
        self.next_append_clears = True
        
    def append_atoms(self, atoms, filename=None):
        "Append an atoms object to the images already stored."
        assert len(atoms) == self.natoms
        if self.next_append_clears:
            i = 0
        else:
            i = self.nimages
        for name in ('P', 'V', 'E', 'K', 'F', 'M', 'A', 'T'):
            a = getattr(self, name)
            newa = np.empty( (i+1,) + a.shape[1:], a.dtype )
            if not self.next_append_clears:
                newa[:-1] = a
            setattr(self, name, newa)
        self.next_append_clears = False
        self.P[i] = atoms.get_positions()
        self.V[i] = atoms.get_velocities()
        self.A[i] = atoms.get_cell()
        try:
            self.E[i] = atoms.get_potential_energy()
        except RuntimeError:
            self.E[i] = np.nan
        self.K[i] = atoms.get_kinetic_energy()
        try:
            self.F[i] = atoms.get_forces(apply_constraint=False)
        except RuntimeError:
            self.F[i] = np.nan
        try:
            self.M[i] = atoms.get_magnetic_moments()
        except (RuntimeError, AttributeError):
            self.M[i] = np.nan
        try:
            self.T[i] = atoms.get_tags()
        except AttributeError:
            if i == 0:
                self.T[i] = 0
            else:
                self.T[i] = self.T[i-1]
        self.nimages = i + 1
        self.filenames.append(filename)
        self.set_dynamic()
        return self.nimages
        
    def set_radii(self, scale):
        if self.shapes == None:
            self.r = self.covalent_radii[self.Z] * scale
        else:
            self.r = np.sqrt(np.sum(self.shapes**2, axis=1)) * scale
                
    def read(self, filenames, index=-1, filetype=None):
        images = []
        names = []
        for filename in filenames:
            i = read(filename, index,filetype)
            
            if not isinstance(i, list):
                i = [i]
            images.extend(i)
            names.extend([filename] * len(i))
            
        self.initialize(images, names)
    
    def import_atoms(self, filename, cur_frame):
        if filename:
            filename = filename[0]
            old_a = self.get_atoms(cur_frame)
            imp_a = read(filename, -1)
            new_a = old_a + imp_a
            self.initialize([new_a], [filename])
    
    def repeat_images(self, repeat):
        n = self.repeat.prod()
        repeat = np.array(repeat)
        self.repeat = repeat
        N = repeat.prod()
        natoms = self.natoms // n
        P = np.empty((self.nimages, natoms * N, 3))
        V = np.empty((self.nimages, natoms * N, 3))
        M = np.empty((self.nimages, natoms * N))
        T = np.empty((self.nimages, natoms * N), int)
        F = np.empty((self.nimages, natoms * N, 3))
        Z = np.empty(natoms * N, int)
        r = np.empty(natoms * N)
        dynamic = np.empty(natoms * N, bool)
        a0 = 0
        for i0 in range(repeat[0]):
            for i1 in range(repeat[1]):
                for i2 in range(repeat[2]):
                    a1 = a0 + natoms
                    for i in range(self.nimages):
                        P[i, a0:a1] = (self.P[i, :natoms] +
                                       np.dot((i0, i1, i2), self.A[i]))
                    V[:, a0:a1] = self.V[:, :natoms]
                    F[:, a0:a1] = self.F[:, :natoms]
                    M[:, a0:a1] = self.M[:, :natoms]
                    T[:, a0:a1] = self.T[:, :natoms]
                    Z[a0:a1] = self.Z[:natoms]
                    r[a0:a1] = self.r[:natoms]
                    dynamic[a0:a1] = self.dynamic[:natoms]
                    a0 = a1
        self.P = P
        self.V = V
        self.F = F
        self.Z = Z
        self.T = T
        self.M = M
        self.r = r
        self.dynamic = dynamic
        self.natoms = natoms * N
        self.selected = np.zeros(natoms * N, bool)
        self.atoms_to_rotate_0 = np.zeros(self.natoms, bool)
        self.visible = np.ones(natoms * N, bool)
        self.nselected = 0

    def center(self):
        """ center each image in the existing unit cell, keeping the cell constant. """
        c = self.A.sum(axis=1) / 2.0 - self.P.mean(axis=1)
        self.P += c[:, np.newaxis, :]
            
    def graph(self, expr):
        """ routine to create the data in ag graphs, defined by the string expr.  """
        import ase.units as units
        code = compile(expr + ',', 'atoms.py', 'eval')

        n = self.nimages
        def d(n1, n2):
            return sqrt(((R[n1] - R[n2])**2).sum())
        def a(n1, n2, n3):
            v1 = R[n1]-R[n2]
            v2 = R[n3]-R[n2]
            arg = np.vdot(v1,v2)/(sqrt((v1**2).sum()*(v2**2).sum()))
            if arg > 1.0: arg = 1.0
            if arg < -1.0: arg = -1.0
            return 180.0*np.arccos(arg)/np.pi
        def dih(n1, n2, n3, n4):
            # vector 0->1, 1->2, 2->3 and their normalized cross products:
            a    = R[n2]-R[n1]
            b    = R[n3]-R[n2]
            c    = R[n4]-R[n3]
            bxa  = np.cross(b,a)
            bxa /= np.sqrt(np.vdot(bxa,bxa))
            cxb  = np.cross(c,b)
            cxb /= np.sqrt(np.vdot(cxb,cxb))
            angle = np.vdot(bxa,cxb)
            # check for numerical trouble due to finite precision:
            if angle < -1: angle = -1
            if angle >  1: angle =  1
            angle = np.arccos(angle)
            if (np.vdot(bxa,c)) > 0: angle = 2*np.pi-angle
            return angle*180.0/np.pi
        # get number of mobile atoms for temperature calculation
        ndynamic = 0
        for dyn in self.dynamic: 
            if dyn: ndynamic += 1
        S = self.selected
        D = self.dynamic[:, np.newaxis]
        E = self.E
        s = 0.0
        data = []
        for i in range(n):
            R = self.P[i]
            V = self.V[i]
            F = self.F[i]
            A = self.A[i]
            M = self.M[i]
            f = ((F * D)**2).sum(1)**.5
            fmax = max(f)
            fave = f.mean()
            epot = E[i]
            ekin = self.K[i]
            e = epot + ekin
            T = 2.0 * ekin / (3.0 * ndynamic * units.kB)
            data = eval(code)
            if i == 0:
                m = len(data)
                xy = np.empty((m, n))
            xy[:, i] = data
            if i + 1 < n:
                s += sqrt(((self.P[i + 1] - R)**2).sum())
        return xy

    def set_dynamic(self, constraints = None):
        self.dynamic = np.ones(self.natoms, bool)
        if constraints is not None:
            for con in constraints: 
                if isinstance(con,FixAtoms):
                    self.dynamic[con.index] = False

    def write(self, filename, rotations='', show_unit_cell=False, bbox=None, **kwargs):
        indices = range(self.nimages)
        p = filename.rfind('@')
        if p != -1:
            try:
                slice = string2index(filename[p + 1:])
            except ValueError:
                pass
            else:
                indices = indices[slice]
                filename = filename[:p]
                if isinstance(indices, int):
                    indices = [indices]

        images = [self.get_atoms(i) for i in indices]
        if len(filename) > 4 and filename[-4:] in ['.eps', '.png', '.pov']:
            write(filename, images, 
                  rotation=rotations, show_unit_cell=show_unit_cell,
                  bbox=bbox, **kwargs)
        else:
            write(filename, images, **kwargs)

    def get_atoms(self, frame):
        atoms = Atoms(positions=self.P[frame],
                      numbers=self.Z,
                      magmoms=self.M[0],
                      tags=self.T[frame],
                      cell=self.A[frame],
                      pbc=self.pbc)

        if not np.isnan(self.V).any():
            atoms.set_velocities(self.V[frame])
        
        # check for constrained atoms and add them accordingly:
        if not self.dynamic.all():
            atoms.set_constraint(FixAtoms(mask=1-self.dynamic))
        
        atoms.set_calculator(SinglePointCalculator(self.E[frame],
                                                   self.F[frame],
                                                   None, None, atoms))
        return atoms
                           
    def delete(self, i):
        self.nimages -= 1
        P = np.empty((self.nimages, self.natoms, 3))
        V = np.empty((self.nimages, self.natoms, 3))
        F = np.empty((self.nimages, self.natoms, 3))
        A = np.empty((self.nimages, 3, 3))
        E = np.empty(self.nimages)
        P[:i] = self.P[:i]
        P[i:] = self.P[i + 1:]
        self.P = P
        V[:i] = self.V[:i]
        V[i:] = self.V[i + 1:]
        self.V = V
        F[:i] = self.F[:i]
        F[i:] = self.F[i + 1:]
        self.F = F
        A[:i] = self.A[:i]
        A[i:] = self.A[i + 1:]
        self.A = A
        E[:i] = self.E[:i]
        E[i:] = self.E[i + 1:]
        self.E = E
        del self.filenames[i]

    def aneb(self):
        n = self.nimages
        assert n % 5 == 0
        levels = n // 5
        n = self.nimages = 2 * levels + 3
        P = np.empty((self.nimages, self.natoms, 3))
        V = np.empty((self.nimages, self.natoms, 3))
        F = np.empty((self.nimages, self.natoms, 3))
        E = np.empty(self.nimages)
        for L in range(levels):
            P[L] = self.P[L * 5]
            P[n - L - 1] = self.P[L * 5 + 4]
            V[L] = self.V[L * 5]
            V[n - L - 1] = self.V[L * 5 + 4]
            F[L] = self.F[L * 5]
            F[n - L - 1] = self.F[L * 5 + 4]
            E[L] = self.E[L * 5]
            E[n - L - 1] = self.E[L * 5 + 4]
        for i in range(3):
            P[levels + i] = self.P[levels * 5 - 4 + i]
            V[levels + i] = self.V[levels * 5 - 4 + i]
            F[levels + i] = self.F[levels * 5 - 4 + i]
            E[levels + i] = self.E[levels * 5 - 4 + i]
        self.P = P
        self.V = V
        self.F = F
        self.E = E

    def interpolate(self, m):
        assert self.nimages == 2
        self.nimages = 2 + m
        P = np.empty((self.nimages, self.natoms, 3))
        V = np.empty((self.nimages, self.natoms, 3))
        F = np.empty((self.nimages, self.natoms, 3))
        A = np.empty((self.nimages, 3, 3))
        E = np.empty(self.nimages)
        P[0] = self.P[0]
        V[0] = self.V[0]
        F[0] = self.F[0]
        A[0] = self.A[0]
        E[0] = self.E[0]
        for i in range(1, m + 1):
            x = i / (m + 1.0)
            y = 1 - x
            P[i] = y * self.P[0] + x * self.P[1]
            V[i] = y * self.V[0] + x * self.V[1]
            F[i] = y * self.F[0] + x * self.F[1]
            A[i] = y * self.A[0] + x * self.A[1]
            E[i] = y * self.E[0] + x * self.E[1]
        P[-1] = self.P[1]
        V[-1] = self.V[1]
        F[-1] = self.F[1]
        A[-1] = self.A[1]
        E[-1] = self.E[1]
        self.P = P
        self.V = V
        self.F = F
        self.A = A
        self.E = E
        self.filenames[1:1] = [None] * m

if __name__ == '__main__':
    import os
    os.system('python gui.py')
