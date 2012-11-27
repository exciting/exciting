import os
import math
import numpy as np
import cPickle as pickle

from ase import Atoms
from ase.data import chemical_symbols
from ase.cluster.base import ClusterBase


class Cluster(Atoms, ClusterBase):
    symmetry = None
    surfaces = None
    lattice_basis = None
    resiproc_basis = None
    atomic_basis = None

    def copy(self):
        cluster = Atoms.copy(self)
        cluster.symmetry = self.symmetry
        cluster.surfaces = self.surfaces.copy()
        cluster.lattice_basis = self.lattice_basis.copy()
        cluster.atomic_basis = self.atomic_basis.copy()
        cluster.resiproc_basis = self.resiproc_basis.copy()
        return cluster

    def get_surfaces(self):
        """Returns the miller indexs of the stored surfaces of the cluster."""
        if not self.surfaces is None:
            return self.surfaces.copy()
        else:
            return None

    def get_layers(self):
        """Return number of atomic layers in stored surfaces directions."""

        layers = []

        for s in self.surfaces:
            n = self.miller_to_direction(s)
            c = self.get_positions().mean(axis=0)
            r = np.dot(self.get_positions() - c, n).max()
            d = self.get_layer_distance(s, 2)
            l = 2 * np.round(r / d).astype(int)

            ls = np.arange(l - 1, l + 2)
            ds = np.array([self.get_layer_distance(s, i) for i in ls])

            mask = (np.abs(ds - r) < 1e-10)

            layers.append(ls[mask][0])

        return np.array(layers, int)

    def get_diameter(self, method='volume'):
        """Returns an estimate of the cluster diameter based on two different
        methods.

        method = 'volume': Returns the diameter of a sphere with the
                           same volume as the atoms. (Default)
        
        method = 'shape': Returns the averaged diameter calculated from the
                          directions given by the defined surfaces.
        """

        if method == 'shape':
            cen = self.get_positions().mean(axis=0)
            pos = self.get_positions() - cen
            d = 0.0
            for s in self.surfaces:
                n = self.miller_to_direction(s)
                r = np.dot(pos, n)
                d += r.max() - r.min()
            return d / len(self.surfaces)
        elif method == 'volume':
            V_cell = np.abs(np.linalg.det(self.lattice_basis))
            N_cell = len(self.atomic_basis)
            N = len(self)
            return 2.0 * (3.0 * N * V_cell /
                          (4.0 * math.pi * N_cell)) ** (1.0 / 3.0)
        else:
            return 0.0

    #Functions to store the cluster
    def write(self, filename=None):
        if not isinstance(filename, str):
            raise Warning('You must specify a valid filename.')

        if os.path.isfile(filename):
            os.rename(filename, filename + '.bak')

        d = {'symmetry': self.symmetry,
             'surfaces': self.surfaces,
             'lattice_basis': self.lattice_basis,
             'resiproc_basis': self.resiproc_basis,
             'atomic_basis': self.atomic_basis,
             'cell': self.get_cell(),
             'pbc': self.get_pbc()}

        f = open(filename, 'w')
        f.write('Cluster')
        pickle.dump(d, f)
        pickle.dump(self.arrays, f)
        f.close()

    def read(self, filename):
        if not os.path.isfile(filename):
            raise Warning('The file specified do not exist.')

        f = open(filename, 'r')

        try:
            if f.read(len('Cluster')) != 'Cluster':
                raise Warning('This is not a compatible file.')
            d = pickle.load(f)
            self.arrays = pickle.load(f)
        except EOFError:
            raise Warinig('Bad file.')

        f.close()

        self.symmetry = d['symmetry']
        self.surfaces = d['surfaces']
        self.lattice_basis = d['lattice_basis']
        self.resiproc_basis = d['resiproc_basis']
        self.atomic_basis = d['atomic_basis']
        self.set_cell(d['cell'])
        self.set_pbc(d['pbc'])
        self.set_constraint()
        self.adsorbate_info = {}
        self.calc = None
