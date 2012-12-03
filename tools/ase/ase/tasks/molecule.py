from math import sqrt
import optparse

import numpy as np

from ase.atoms import Atoms, string2symbols
from ase.structure import molecule
from ase.tasks.task import OptimizeTask
from ase.data import covalent_radii, atomic_numbers
from ase.data import ground_state_magnetic_moments
from ase.utils.eos import EquationOfState
from ase.io.trajectory import PickleTrajectory
import ase.units as units


class MoleculeTask(OptimizeTask):
    taskname = 'molecule'

    def __init__(self, vacuum=3.0, cell=None, atomize=False,
                 bond_length=None, fit=None,
                 **kwargs):
        """Molecule task.

        This task can calculate bond lengths and vibration frequencies
        of dimer molecules."""

        self.vacuum = vacuum
        self.unit_cell = cell
        self.atomize = atomize
        self.bond_length = bond_length
        self.fit = fit

        OptimizeTask.__init__(self, **kwargs)

        self.summary_keys = ['energy', 'relaxed energy', 'distance',
                             'frequency', 'atomic energy']

    def run(self, names1):
        names = []
        atoms = set()
        for name in names1:
            if name.lower() == 'g2':
                from ase.data.g2 import molecule_names
                names.extend(molecule_names)
                from ase.data.g2 import atom_names
                if self.atomize:
                    atoms.update(atom_names)
            elif name.lower() == 'g2-1':
                from ase.data.g2_1 import molecule_names
                names.extend(molecule_names)
                from ase.data.g2_1 import atom_names
                if self.atomize:
                    atoms.update(atom_names)
            else:
                names.append(name)
                if self.atomize:
                    atoms.update(self.build_system(name).get_chemical_symbols())
        if self.atomize:
            names.extend(atoms)

        return OptimizeTask.run(self, names)

    def build_system(self, name):
        try:
            # Known molecule or atom?
            atoms = molecule(name)
            if len(atoms) == 2 and self.bond_length is not None:
                atoms.set_distance(0, 1, self.bond_length)
        except NotImplementedError:
            symbols = string2symbols(name)
            if len(symbols) == 1:
                Z = atomic_numbers[symbols[0]]
                magmom = ground_state_magnetic_moments[Z]
                atoms = Atoms(name, magmoms=[magmom])
            elif len(symbols) == 2:
                # Dimer
                if self.bond_length is None:
                    b = (covalent_radii[atomic_numbers[symbols[0]]] +
                         covalent_radii[atomic_numbers[symbols[1]]])
                else:
                    b = self.bond_length
                atoms = Atoms(name, positions=[(0, 0, 0),
                                               (b, 0, 0)])
            else:
                raise ValueError('Unknown molecule: ' + name)

        if self.unit_cell is None:
            atoms.center(vacuum=self.vacuum)
        else:
            atoms.cell = self.unit_cell
            atoms.center()

        return atoms

    def fit_bond_length(self, name, atoms, data):
        N, x = self.fit
        assert N % 2 == 1
        d0 = atoms.get_distance(0, 1)
        distances = np.linspace(d0 * (1 - x), d0 * (1 + x), N)
        energies = []
        traj = PickleTrajectory(self.get_filename(name, 'fit.traj'), 'w')
        for d in distances:
            atoms.set_distance(0, 1, d)
            energies.append(atoms.get_potential_energy())
            self.check_occupation_numbers(atoms)
            traj.write(atoms)

        traj.close()

        data['distances'] = distances
        data['energies'] = energies

    def calculate(self, name, atoms):
        data = OptimizeTask.calculate(self, name, atoms)
        if self.fmax is not None and len(atoms) == 2:
            data['distance'] = atoms.get_distance(0, 1)
        self.check_occupation_numbers(atoms)
        if self.fit and len(atoms) == 2:
            self.fit_bond_length(name, atoms, data)
        return data

    def analyse(self):
        for name, data in self.data.items():
            if 'distances' in data:
                distances = data['distances']
                energies = data['energies']
                fit0 = np.poly1d(np.polyfit(1 / distances, energies, 3))
                fit1 = np.polyder(fit0, 1)
                fit2 = np.polyder(fit1, 1)

                dmin = None
                for t in np.roots(fit1):
                    if t > 0 and fit2(t) > 0:
                        dmin = 1 / t
                        break

                if dmin is None:
                    raise ValueError('No minimum!')

                if abs(dmin) < min(distances) or abs(dmin) > max(distances):
                    raise ValueError(name + ': fit outside of range! ' + \
                                     str(abs(dmin)) + ' not in ' + \
                                     str(distances))

                emin = fit0(t)
                k = fit2(t) * t**4
                m1, m2 = self.create_system(name).get_masses()
                m = m1 * m2 / (m1 + m2)
                hnu = units._hbar * 1e10 * sqrt(k / units._e / units._amu / m)

                data['relaxed energy'] = emin
                data['distance'] = dmin
                data['frequency'] = hnu

        for name, data in self.data.items():
            atoms = self.create_system(name)
            if len(atoms) == 1:
                continue

            eatoms = 0.0
            for symbol in atoms.get_chemical_symbols():
                if symbol in self.data and symbol != name:
                    eatoms += self.data[symbol]['energy']
                else:
                    eatoms = None
                    break
            ea = None
            ea0 = None
            if eatoms is not None:
                data['atomic energy'] = eatoms

    def add_options(self, parser):
        OptimizeTask.add_options(self, parser)

        mol = optparse.OptionGroup(parser, 'Molecule')
        mol.add_option('-v', '--vacuum', type='float', default=3.0,
                       help='Amount of vacuum to add around isolated systems '
                       '(in Angstrom).')
        mol.add_option('--unit-cell',
                       help='Unit cell.  Examples: "10.0" or "9,10,11" ' +
                       '(in Angstrom).')
        mol.add_option('--bond-length', type='float',
                       help='Bond length of dimer in Angstrom.')
        mol.add_option('-F', '--fit', metavar='N,x',
                       help='Find optimal bondlength and vibration ' +
                       'frequency using N points and displacements from ' +
                       '-x % to +x %.')
        mol.add_option('--atomize', action='store_true',
                       help='Calculate Atomization energies.')
        parser.add_option_group(mol)

    def parse(self, opts, args):
        OptimizeTask.parse(self, opts, args)

        self.vacuum = opts.vacuum
        self.bond_length = opts.bond_length
        self.atomize = opts.atomize

        if opts.fit:
            points, strain = opts.fit.split(',')
            self.fit = (int(points), float(strain) * 0.01)

        if opts.unit_cell:
            if ',' in opts.unit_cell:
                self.unit_cell = [float(x) for x in opts.unit_cell.split(',')]
            else:
                self.unit_cell = [float(opts.unit_cell)] * 3

    def check_occupation_numbers(self, config):
        """Check that occupation numbers are integers and sum
        to desired magnetic moment.  """
        if config.pbc.any():
            return
        calc = config.get_calculator()
        try:
            mref = abs(config.get_initial_magnetic_moments().sum())
            nspins = calc.get_number_of_spins()
            mcalc = 0.0
            for s in range(nspins):
                f = calc.get_occupation_numbers(spin=s)
                if abs((f.round() - f).sum()) > 0.0001:
                    raise RuntimeError('Fractional occupation numbers?! ' + \
                                       str(f) + ' for spin ' + str(s))
                mcalc += abs(f).sum() * (-1)**s
            mcalc = abs(mcalc)
            if mref > 0.0:
                if  abs(mcalc - mref) > 0.0001:
                    raise RuntimeError('Incorrect magnetic moment?! ' + \
                                       str(mcalc) + ' vs ' + str(mref))

        except AttributeError:
            pass
