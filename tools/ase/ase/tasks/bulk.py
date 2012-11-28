import optparse

import numpy as np

from ase.lattice import bulk
from ase.tasks.task import OptimizeTask
from ase.data import chemical_symbols, reference_states
from ase.utils.eos import EquationOfState
from ase.io.trajectory import PickleTrajectory


class BulkTask(OptimizeTask):
    taskname = 'bulk'

    def __init__(self, crystal_structure=None, lattice_constant=None,
                 c_over_a=None, cubic=False, orthorhombic=False, fit=None,
                 **kwargs):
        """Bulk task."""

        self.crystal_structure = crystal_structure
        self.lattice_constant = lattice_constant
        self.c_over_a = c_over_a
        self.cubic = cubic
        self.orthorhombic = orthorhombic
        self.fit = fit

        self.repeat = None

        OptimizeTask.__init__(self, **kwargs)

        self.summary_keys = ['energy', 'fitted energy', 'volume', 'B']

    def expand(self, names):
        """Expand fcc, bcc, hcp and diamond.

        The name fcc will be expanded to all the elements with the fcc
        stucture and so on."""

        names = OptimizeTask.expand(self, names)

        newnames = []
        for name in names:
            if name in ['fcc', 'bcc', 'hcp', 'diamond']:
                for Z in range(1, 95):
                    x = reference_states[Z]
                    if x is not None and x['symmetry'] == name:
                        newnames.append(chemical_symbols[Z])
            else:
                newnames.append(name)

        return newnames

    def build_system(self, name):
        atoms = bulk(name, crystalstructure=self.crystal_structure,
                     a=self.lattice_constant, covera=self.c_over_a,
                     orthorhombic=self.orthorhombic, cubic=self.cubic)

        M = {'Fe': 2.3, 'Co': 1.2, 'Ni': 0.6}.get(name)
        if M is not None:
            atoms.set_initial_magnetic_moments([M] * len(atoms))

        if self.repeat is not None:
            r = self.repeat.split(',')
            if len(r) == 1:
                r = 3 * r
            atoms = atoms.repeat([int(c) for c in r])

        return atoms

    def fit_volume(self, name, atoms):
        N, x = self.fit
        cell0 = atoms.get_cell()
        strains = np.linspace(1 - x, 1 + x, N)
        energies = []
        traj = PickleTrajectory(self.get_filename(name, 'fit.traj'), 'w')
        for s in strains:
            atoms.set_cell(cell0 * s, scale_atoms=True)
            energies.append(atoms.get_potential_energy())
            traj.write(atoms)

        traj.close()

        assert N % 2 == 1
        data = {'energy': energies[N // 2],
                'strains': strains,
                'energies': energies}

        return data

    def calculate(self, name, atoms):
        #????
        if self.fit:
            return self.fit_volume(name, atoms)
        else:
            return OptimizeTask.calculate(self, name, atoms)

    def analyse(self):
        for name, data in self.data.items():
            if 'strains' in data:
                atoms = self.create_system(name)
                volumes = data['strains']**3 * atoms.get_volume()
                energies = data['energies']
                eos = EquationOfState(volumes, energies)
                try:
                    v, e, B = eos.fit()
                except ValueError:
                    pass
                else:
                    data['fitted energy'] = e
                    data['volume'] = v
                    data['B'] = B

                    if abs(v) < min(volumes) or abs(v) > max(volumes):
                        raise ValueError(name + ': fit outside of range! ' + \
                                         str(abs(v)) + ' not in ' + \
                                         str(volumes))

    def add_options(self, parser):
        OptimizeTask.add_options(self, parser)

        bulk = optparse.OptionGroup(parser, 'Bulk')
        bulk.add_option('-F', '--fit', metavar='N,x',
                        help='Find optimal volume and bulk modulus ' +
                        'using N points and variations of the lattice ' +
                        'constants from -x % to +x %.')
        bulk.add_option('-x', '--crystal-structure',
                        help='Crystal structure.',
                        choices=['sc', 'fcc', 'bcc', 'diamond', 'hcp',
                                 'zincblende', 'rocksalt',
                                 'cesiumchloride', 'fluorite'])
        bulk.add_option('-a', '--lattice-constant', type='float',
                        help='Lattice constant in Angstrom.')
        bulk.add_option('--c-over-a', type='float',
                        help='c/a ratio.')
        bulk.add_option('-O', '--orthorhombic', action='store_true',
                        help='Use orthorhombic unit cell.')
        bulk.add_option('-C', '--cubic', action='store_true',
                        help='Use cubic unit cell.')
        bulk.add_option('-r', '--repeat',
                        help='Repeat unit cell.  Use "-r 2" or "-r 2,3,1".')
        parser.add_option_group(bulk)

    def parse(self, opts, args):
        OptimizeTask.parse(self, opts, args)

        if opts.fit:
            points, strain = opts.fit.split(',')
            self.fit = (int(points), float(strain) * 0.01)

        self.crystal_structure = opts.crystal_structure
        self.lattice_constant = opts.lattice_constant
        self.c_over_a = opts.c_over_a
        self.orthorhombic = opts.orthorhombic
        self.cubic = opts.cubic
        self.repeat = opts.repeat
