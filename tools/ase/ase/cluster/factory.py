import numpy as np

from ase.data import atomic_numbers as ref_atomic_numbers
from ase.lattice.spacegroup import Spacegroup
from ase.cluster.base import ClusterBase
from ase.cluster.cluster import Cluster

class ClusterFactory(ClusterBase):
    directions = [[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]]

    atomic_basis = np.array([[0., 0., 0.]])

    element_basis = None

    def __call__(self, symbols, surfaces, layers, latticeconstant=None,
                 center=None, vacuum=0.0, debug=0):
        self.debug = debug

        # Interpret symbol
        self.set_atomic_numbers(symbols)

        # Interpret lattice constant
        if latticeconstant is None:
            if self.element_basis is None:
                self.lattice_constant = self.get_lattice_constant()
            else:
                raise ValueError("A lattice constant must be specified for a compound")
        else:
            self.lattice_constant = latticeconstant

        self.set_basis()

        if self.debug:
            print "Lattice constant(s):", self.lattice_constant
            print "Lattice basis:\n", self.lattice_basis
            print "Resiprocal basis:\n", self.resiproc_basis
            print "Atomic basis:\n", self.atomic_basis

        self.set_surfaces_layers(surfaces, layers)
        self.set_lattice_size(center)

        if self.debug:
            print "Center position:", self.center.round(2)
            print "Base lattice size:", self.size

        cluster = self.make_cluster(vacuum)
        cluster.symmetry = self.xtal_name
        cluster.surfaces = self.surfaces.copy()
        cluster.lattice_basis = self.lattice_basis.copy()
        cluster.atomic_basis = self.atomic_basis.copy()
        cluster.resiproc_basis = self.resiproc_basis.copy()
        return cluster

    def make_cluster(self, vacuum):
        # Make the base crystal by repeating the unit cell
        size = np.array(self.size)
        translations = np.zeros((size.prod(), 3))
        for h in range(size[0]):
            for k in range(size[1]):
                for l in range(size[2]):
                    i = h * (size[1] * size[2]) + k * size[2] + l
                    translations[i] = np.dot([h, k, l], self.lattice_basis)

        atomic_basis = np.dot(self.atomic_basis, self.lattice_basis)
        positions = np.zeros((len(translations) * len(atomic_basis), 3))
        numbers = np.zeros(len(positions))
        n = len(atomic_basis)
        for i, trans in enumerate(translations):
            positions[n*i:n*(i+1)] = atomic_basis + trans
            numbers[n*i:n*(i+1)] = self.atomic_numbers

        # Remove all atoms that is outside the defined surfaces
        for s, l in zip(self.surfaces, self.layers):
            n = self.miller_to_direction(s)
            rmax = self.get_layer_distance(s, l + 0.1)

            r = np.dot(positions - self.center, n)
            mask = np.less(r, rmax)

            if self.debug > 1:
                print "Cutting %s at %i layers ~ %.3f A" % (s, l, rmax)

            positions = positions[mask]
            numbers = numbers[mask]

        # Fit the cell, so it only just consist the atoms
        min = np.zeros(3)
        max = np.zeros(3)
        for i in range(3):
            v = self.directions[i]
            r = np.dot(positions, v)
            min[i] = r.min()
            max[i] = r.max()

        cell = max - min + vacuum
        positions = positions - min + vacuum / 2.0
        self.center = self.center - min + vacuum / 2.0

        return Cluster(symbols=numbers, positions=positions, cell=cell)

    def set_atomic_numbers(self, symbols):
        "Extract atomic number from element"
        # The types that can be elements: integers and strings
        atomic_numbers = []
        if self.element_basis is None:
            if isinstance(symbols, str):
                atomic_numbers.append(ref_atomic_numbers[symbols])
            elif isinstance(symbols, int):
                atomic_numbers.append(symbols)
            else:
                raise TypeError("The symbol argument must be a " +
                                "string or an atomic number.")
            element_basis = [0] * len(self.atomic_basis)
        else:
            if isinstance(symbols, (list, tuple)):
                nsymbols = len(symbols)
            else:
                nsymbols = 0

            nelement_basis = max(self.element_basis) + 1
            if nsymbols != nelement_basis:
                raise TypeError("The symbol argument must be a sequence " +
                                "of length %d" % (nelement_basis,) +
                                " (one for each kind of lattice position")

            for s in symbols:
                if isinstance(s, str):
                    atomic_numbers.append(ref_atomic_numbers[s])
                elif isinstance(s, int):
                    atomic_numbers.append(s)
                else:
                    raise TypeError("The symbol argument must be a " +
                                    "string or an atomic number.")
            element_basis = self.element_basis

        self.atomic_numbers = [atomic_numbers[n] for n in element_basis]
        assert len(self.atomic_numbers) == len(self.atomic_basis)

    def set_lattice_size(self, center):
        if center is None:
            offset = np.zeros(3)
        else:
            offset = np.array(center)
            if (offset > 1.0).any() or (offset < 0.0).any():
                raise ValueError("Center offset must lie within the lattice unit \
                                  cell.")

        max = np.ones(3)
        min = -np.ones(3)
        v = np.linalg.inv(self.lattice_basis.T)
        for s, l in zip(self.surfaces, self.layers):
            n = self.miller_to_direction(s) * self.get_layer_distance(s, l)
            k = np.round(np.dot(v, n), 2)
            for i in range(3):
                if k[i] > 0.0:
                    k[i] = np.ceil(k[i])
                elif k[i] < 0.0:
                    k[i] = np.floor(k[i])

            if self.debug > 1:
                print "Spaning %i layers in %s in lattice basis ~ %s" % (l, s, k)

            max[k > max] = k[k > max]
            min[k < min] = k[k < min]

        self.center = np.dot(offset - min, self.lattice_basis)
        self.size = (max - min + np.ones(3)).astype(int)

    def set_surfaces_layers(self, surfaces, layers):
        if len(surfaces) != len(layers):
            raise ValueError("Improper size of surface and layer arrays: %i != %i"
                             % (len(surfaces), len(layers)))

        sg = Spacegroup(self.spacegroup)
        surfaces = np.array(surfaces)
        layers = np.array(layers)

        for i, s in enumerate(surfaces):
            s = reduce_miller(s)
            surfaces[i] = s

        surfaces_full = surfaces.copy()
        layers_full = layers.copy()

        for s, l in zip(surfaces, layers):
            equivalent_surfaces = sg.equivalent_reflections(s.reshape(-1, 3))

            for es in equivalent_surfaces:
                # If the equivalent surface (es) is not in the surface list,
                # then append it.
                if not np.equal(es, surfaces_full).all(axis=1).any():
                    surfaces_full = np.append(surfaces_full, es.reshape(1, 3), axis=0)
                    layers_full = np.append(layers_full, l)

        self.surfaces = surfaces_full.copy()
        self.layers = layers_full.copy()

    def get_resiproc_basis(self, basis):
        """Returns the resiprocal basis to a given lattice (crystal) basis"""
        k = 1 / np.dot(basis[0], cross(basis[1], basis[2]))

        # The same as the inversed basis matrix transposed
        return k * np.array([cross(basis[1], basis[2]),
                             cross(basis[2], basis[0]),
                             cross(basis[0], basis[1])])

# Helping functions
def cross(a, b):
    """The cross product of two vectors."""
    return np.array([a[1]*b[2] - b[1]*a[2],
                     a[2]*b[0] - b[2]*a[0],
                     a[0]*b[1] - b[0]*a[1]])

def GCD(a,b):
    """Greatest Common Divisor of a and b."""
    #print "--"
    while a != 0:
        #print a,b,">",
        a,b = b%a,a
        #print a,b
    return b

def reduce_miller(hkl):
    """Reduce Miller index to the lowest equivalent integers."""
    hkl = np.array(hkl)
    old = hkl.copy()

    d = GCD(GCD(hkl[0], hkl[1]), hkl[2])
    while d != 1:
        hkl = hkl / d
        d = GCD(GCD(hkl[0], hkl[1]), hkl[2])

    if np.dot(old, hkl) > 0:
        return hkl
    else:
        return -hkl

