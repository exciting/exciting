from math import sqrt

import numpy as np

import ase.parallel as mpi
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read


class NEB:
    def __init__(self, images, k=0.1, climb=False, parallel=False,
                 world=None):
        """Nudged elastic band.

        images: list of Atoms objects
            Images defining path from initial to final state.
        k: float or list of floats
            Spring constant(s).  One number or one for each spring.
        climb: bool
            Use a climbing image (default is no climbing image).
        parallel: bool
            Distribute images over processors.
        """
        self.images = images
        self.climb = climb
        self.parallel = parallel
        self.natoms = len(images[0])
        self.nimages = len(images)
        self.emax = np.nan

        if isinstance(k, (float, int)):
            k = [k] * (self.nimages - 1)
        self.k = list(k)

        if world is None:
            world = mpi.world
        self.world = world

        assert not parallel or world.size % (self.nimages - 2) == 0

    def interpolate(self):
        pos1 = self.images[0].get_positions()
        pos2 = self.images[-1].get_positions()
        d = (pos2 - pos1) / (self.nimages - 1.0)
        for i in range(1, self.nimages - 1):
            self.images[i].set_positions(pos1 + i * d)
            # Parallel NEB with Jacapo needs this:
            try:
                self.images[i].get_calculator().set_atoms(self.images[i])
            except AttributeError:
                pass
            
    def get_positions(self):
        positions = np.empty(((self.nimages - 2) * self.natoms, 3))
        n1 = 0
        for image in self.images[1:-1]:
            n2 = n1 + self.natoms
            positions[n1:n2] = image.get_positions()
            n1 = n2
        return positions

    def set_positions(self, positions):
        n1 = 0
        for image in self.images[1:-1]:
            n2 = n1 + self.natoms
            image.set_positions(positions[n1:n2])
            n1 = n2

            # Parallel NEB with Jacapo needs this:
            try:
                image.get_calculator().set_atoms(image)
            except AttributeError:
                pass
        
    def get_forces(self):
        """Evaluate and return the forces."""
        images = self.images
        forces = np.empty(((self.nimages - 2), self.natoms, 3))
        energies = np.empty(self.nimages - 2)

        if not self.parallel:
            # Do all images - one at a time:
            for i in range(1, self.nimages - 1):
                energies[i - 1] = images[i].get_potential_energy()
                forces[i - 1] = images[i].get_forces()
        else:
            # Parallelize over images:
            i = self.world.rank * (self.nimages - 2) // self.world.size + 1
            try:
                energies[i - 1] = images[i].get_potential_energy()
                forces[i - 1] = images[i].get_forces()
            except:
                # Make sure other images also fail:
                error = self.world.sum(1.0)
                raise
            else:
                error = self.world.sum(0.0)
                if error:
                    raise RuntimeError('Parallel NEB failed!')
                
            for i in range(1, self.nimages - 1):
                root = (i - 1) * self.world.size // (self.nimages - 2)
                self.world.broadcast(energies[i - 1:i], root)
                self.world.broadcast(forces[i - 1], root)

        imax = 1 + np.argsort(energies)[-1]
        self.emax = energies[imax - 1]
        
        tangent1 = images[1].get_positions() - images[0].get_positions()
        for i in range(1, self.nimages - 1):
            tangent2 = (images[i + 1].get_positions() -
                        images[i].get_positions())
            if i < imax:
                tangent = tangent2
            elif i > imax:
                tangent = tangent1
            else:
                tangent = tangent1 + tangent2
                
            tt = np.vdot(tangent, tangent)
            f = forces[i - 1]
            ft = np.vdot(f, tangent)
            if i == imax and self.climb:
                f -= 2 * ft / tt * tangent
            else:
                f -= ft / tt * tangent
                f -= np.vdot(tangent1 * self.k[i - 1] -
                             tangent2 * self.k[i], tangent) / tt * tangent
                
            tangent1 = tangent2

        return forces.reshape((-1, 3))

    def get_potential_energy(self):
        return self.emax

    def __len__(self):
        return (self.nimages - 2) * self.natoms


class SingleCalculatorNEB(NEB):
    def __init__(self, images, k=0.1, climb=False):
        if isinstance(images, str):
            # this is a filename
            traj = read(images, '0:')
            images = []
            for atoms in traj:
                images.append(atoms)

        NEB.__init__(self, images, k, climb, False)
        self.calculators = [None] * self.nimages
        self.energies_ok = False
 
    def interpolate(self, initial=0, final=-1):
        """Interpolate linearly between initial and final images."""
        if final < 0:
            final = self.nimages + final
        n = final - initial
        pos1 = self.images[initial].get_positions()
        pos2 = self.images[final].get_positions()
        d = (pos2 - pos1) / n
        for i in range(1, n):
            self.images[initial + i].set_positions(pos1 + i * d)

    def refine(self, steps=1, begin=0, end=-1):
        """Refine the NEB trajectory."""
        if end < 0:
            end = self.nimages + end
        j = begin
        n = end - begin
        for i in range(n):
            for k in range(steps):
                self.images.insert(j + 1, self.images[j].copy())
                self.calculators.insert(j + 1, None)
            self.k[j:j + 1] = [self.k[j] * (steps + 1)] * (steps + 1)
            self.nimages = len(self.images)
            self.interpolate(j, j + steps + 1)
            j += steps + 1

    def set_positions(self, positions):
        # new positions -> new forces
        if self.energies_ok:
            # restore calculators
            self.set_calculators(self.calculators[1:-1])
        NEB.set_positions(self, positions)

    def get_calculators(self):
        """Return the original calculators."""
        calculators = []
        for i, image in enumerate(self.images):
            if self.calculators[i] is None:
                calculators.append(image.get_calculator())
            else:
                calculators.append(self.calculators[i])
        return calculators
    
    def set_calculators(self, calculators):
        """Set new calculators to the images."""
        self.energies_ok = False

        if not isinstance(calculators, list):
            calculators = [calculators] * self.nimages

        n = len(calculators)
        if n == self.nimages:
            for i in range(self.nimages):
                self.images[i].set_calculator(calculators[i])
        elif n == self.nimages - 2:
            for i in range(1, self.nimages - 1):
                self.images[i].set_calculator(calculators[i - 1])
        else:
            raise RuntimeError(
                'len(calculators)=%d does not fit to len(images)=%d'
                % (n, self.nimages))

    def get_energies_and_forces(self, all=False):
        """Evaluate energies and forces and hide the calculators"""
        if self.energies_ok:
            return

        self.emax = -1.e32

        def calculate_and_hide(i):
            image = self.images[i]
            calc = image.get_calculator()
            if self.calculators[i] is None:
                self.calculators[i] = calc
            if calc is not None:
                if not isinstance(calc, SinglePointCalculator):
                    self.images[i].set_calculator(
                        SinglePointCalculator(image.get_potential_energy(),
                                              image.get_forces(),
                                              None,
                                              None,
                                              image))
                self.emax = min(self.emax, image.get_potential_energy())

        if all and self.calculators[0] is None:
            calculate_and_hide(0)

        # Do all images - one at a time:
        for i in range(1, self.nimages - 1):
            calculate_and_hide(i)

        if all and self.calculators[-1] is None:
            calculate_and_hide(-1)

        self.energies_ok = True
       
    def get_forces(self):
        self.get_energies_and_forces()
        return NEB.get_forces(self)

    def n(self):
        return self.nimages

    def write(self, filename):
        from ase.io.trajectory import PickleTrajectory
        traj = PickleTrajectory(filename, 'w', self)
        traj.write()
        traj.close()

    def __add__(self, other):
        for image in other:
            self.images.append(image)
        return self


def fit(images):
    E = [i.get_potential_energy() for i in images]
    F = [i.get_forces() for i in images]
    R = [i.get_positions() for i in images]
    return fit0(E, F, R)


def fit0(E, F, R):
    E = np.array(E) - E[0]
    n = len(E)
    Efit = np.empty((n - 1) * 20 + 1)
    Sfit = np.empty((n - 1) * 20 + 1)

    s = [0]
    for i in range(n - 1):
        s.append(s[-1] + sqrt(((R[i + 1] - R[i])**2).sum()))

    lines = []
    for i in range(n):
        if i == 0:
            d = R[1] - R[0]
            ds = 0.5 * s[1]
        elif i == n - 1:
            d = R[-1] - R[-2]
            ds = 0.5 * (s[-1] - s[-2])
        else:
            d = R[i + 1] - R[i - 1]
            ds = 0.25 * (s[i + 1] - s[i - 1])

        d = d / sqrt((d**2).sum())
        dEds = -(F[i] * d).sum()
        x = np.linspace(s[i] - ds, s[i] + ds, 3)
        y = E[i] + dEds * (x - s[i])
        lines.append((x, y))

        if i > 0:
            s0 = s[i - 1]
            s1 = s[i]
            x = np.linspace(s0, s1, 20, endpoint=False)
            c = np.linalg.solve(np.array([(1, s0,   s0**2,     s0**3),
                                          (1, s1,   s1**2,     s1**3),
                                          (0,  1,  2 * s0, 3 * s0**2),
                                          (0,  1,  2 * s1, 3 * s1**2)]),
                                np.array([E[i - 1], E[i], dEds0, dEds]))
            y = c[0] + x * (c[1] + x * (c[2] + x * c[3]))
            Sfit[(i - 1) * 20:i * 20] = x
            Efit[(i - 1) * 20:i * 20] = y
        
        dEds0 = dEds

    Sfit[-1] = s[-1]
    Efit[-1] = E[-1]
    return s, E, Sfit, Efit, lines
