from math import exp, sqrt

import numpy as np

from ase.atoms import Atoms


class STM:
    def __init__(self, atoms, symmetries=None):
        if isinstance(atoms, Atoms):
            calc = atoms.get_calculator()
        else:
            calc = atoms
            atoms = calc.get_atoms()
        self.nbands = calc.get_number_of_bands()
        self.weights = calc.get_k_point_weights()
        self.nkpts = len(self.weights)
        self.nspins = calc.get_number_of_spins()
        self.eigs = np.array([[calc.get_eigenvalues(k, s)
                               for k in range(self.nkpts)]
                              for s in range(self.nspins)])
        self.eigs -= calc.get_fermi_level()
        self.calc = calc
        self.cell = atoms.get_cell()
        assert not self.cell[2, :2].any() and not self.cell[:2, 2].any()
        self.ldos = None
        self.symmetries = symmetries or []
                               
    def calculate_ldos(self, width=None):
        if self.ldos is not None and width == self.width:
            return

        if width is None:
            width = 0.1
            
        ldos = None
        for s in range(self.nspins):
            for k in range(self.nkpts):
                for n in range(self.nbands):
                    psi = self.calc.get_pseudo_wave_function(n, k, s)
                    if ldos is None:
                        ldos = np.zeros(psi.shape)
                    f = (exp(-(self.eigs[s, k, n] / width)**2) *
                         self.weights[k])
                    ldos += f * (psi * np.conj(psi)).real

        if 0 in self.symmetries:
            # (x,y) -> (-x,y)
            ldos[1:] += ldos[:0:-1].copy()
            ldos[1:] *= 0.5

        if 1 in self.symmetries:
            # (x,y) -> (x,-y)
            ldos[:, 1:] += ldos[:, :0:-1].copy()
            ldos[:, 1:] *= 0.5
            
        if 2 in self.symmetries:
            # (x,y) -> (y,x)
            ldos += ldos.transpose((1, 0, 2)).copy()
            ldos *= 0.5
            
        self.ldos = ldos
        self.width = width

    #def save_ldos(self, filename='ldos.pckl'):
        

    def get_averaged_current(self, z, width=None):
        self.calculate_ldos(width)
        nz = self.ldos.shape[2]

        # Find grid point:
        n = z / self.cell[2, 2] * nz
        dn = n - np.floor(n)
        n = int(n) % nz
        print n,dn

        # Average and do linear interpolation:
        return ((1 - dn) * self.ldos[:, :, n].mean() +
                dn *       self.ldos[:, :, (n + 1) % nz].mean())
    
    def scan(self, current, z=None, width=None):
        self.calculate_ldos(width)

        L = self.cell[2, 2]
        if z is None:
            z = L / 2

        nz = self.ldos.shape[2]
        n = int(round(z / L * nz)) % nz
        h = L / nz

        ldos = self.ldos.reshape((-1, nz))

        heights = np.empty(ldos.shape[0])
        for i, a in enumerate(ldos):
            heights[i], z, n = find_height(a, current, z, n, nz, h)

        heights.shape = self.ldos.shape[:2]
        return heights
    
    def linescan(self, current, p1, p2, npoints=None, z=None, width=None):
        self.calculate_ldos(width)

        L = self.cell[2, 2]
        if z is None:
            z = L / 2

        nz = self.ldos.shape[2]
        n = int(round(z / L * nz)) % nz
        h = L / nz
        ldos = self.ldos.reshape((-1, nz))

        p1 = np.asarray(p1)
        p2 = np.asarray(p2)
        d = p2 - p1
        s = sqrt(np.dot(d, d))
        
        if npints == None:
            npoints = int(3 * s / h + 2)

        cell = self.cell[:2, :2]
        shape = np.array(self.ldos.shape[:2], float)
        M = cell.I
        heights = np.empty(npoints)
        for i in range(npoints):
            p = p1 + i * d / (npoints - 1)
            q = np.dot(M, p) * shape
            qi = q.astype(int)
            n0, n1 = qi
            f = q - qi
            g = 1 - f
            a = (g[0] * g[0] * ldos[n0,     n1    ] +
                 f[0] * g[0] * ldos[n0 + 1, n1    ] +
                 g[0] * f[0] * ldos[n0,     n1 + 1] +
                 f[0] * f[0] * ldos[n0 + 1, n1 + 1])
            heights[i], z, n = find_height(a, current, z, n, nz, h)
        return np.linspace(0, s, npoints), heights

    def cube(self, filename, atoms=None):
        pass


def find_height(array, current, z, n, nz, h):
    c1 = array[n]
    sign = cmp(c1, current)
    m = 0
    while m < nz:
        n = (n + sign) % nz
        z += sign * h
        c2 = array[n]
        if cmp(c2, current) != sign:
            break
        c1 = c2
        m += 1

    if m == nz:
        print z, n, nz, h, current, array
        raise RuntimeError('Tip crash!')

    return z - sign * h * (current - c2) / (c1 - c2), z, n

                
            
        
    
        
