#!/usr/bin/env python

"""This module defines extensible classes for running tests on molecules.

Use this to compare different calculators, XC functionals and so on by
calculating e.g. atomization energies and bond lengths across the
g2 database of small molecules.
"""

import os
import sys
import traceback

import numpy as np

from ase import PickleTrajectory, read
from ase.calculators.emt import EMT
from ase.data.molecules import molecule, atoms as g2_atoms, g1


class BatchTest:
    """Contains logic for looping over tests and file management."""
    def __init__(self, test):
        self.test = test
        self.txt = sys.stdout # ?

    def run_single_test(self, formula):
        print >> self.txt, self.test.name, formula, '...',
        self.txt.flush()
        filename = self.test.get_filename(formula)
        if os.path.exists(filename):
            print >> self.txt, 'Skipped.'
            return
        try:
            open(filename, 'w').close() # Empty file
            system, calc = self.test.setup(formula)
            self.test.run(formula, system, filename)
            print >> self.txt, 'OK!'
            self.txt.flush()
        except self.test.exceptions:
            print >> self.txt, 'Failed!'
            traceback.print_exc(file=self.txt)
            print >> self.txt
            self.txt.flush()

    def run(self, formulas):
        """Run a batch of tests.

        This will invoke the test method on each formula, printing
        status to stdout.

        Those formulas that already have test result files will
        be skipped."""
        
        # Create directories if necessary
        if self.test.dir and not os.path.isdir(self.test.dir):
            os.mkdir(self.test.dir) # Won't work on 'dir1/dir2', but oh well
        
        for formula in formulas:
            self.run_single_test(formula)
    
    def collect(self, formulas, verbose=False):
        """Yield results of previous calculations."""
        for formula in formulas:
            try:
                filename = self.test.get_filename(formula)
                results = self.test.retrieve(formula, filename)
                if verbose:
                    print >> self.txt, 'Loaded:', formula, filename
                yield formula, results
            except (IOError, RuntimeError, TypeError):
                # XXX which errors should we actually catch?
                if verbose:
                    print >> self.txt, 'Error:', formula, '[%s]' % filename
                    traceback.print_exc(file=self.txt)


class MoleculeTest:
    """Generic class for runnings various tests on the g2 dataset.

    Usage: instantiate MoleculeTest with desired test settings and
    invoke its run() method on the desired formulas.

    This class will use the ASE EMT calculator by default.  You can
    create a subclass using an arbitrary calculator by overriding the
    setup_calculator method.  Most methods can be overridden to
    provide highly customized behaviour.  """
    
    def __init__(self, name, vacuum=6.0, exceptions=None):
        """Create a molecule test.

        The name parameter will be part of all output files generated
        by this molecule test.  If name contains a '/' character, the
        preceding part will be interpreted as a directory in which to
        put files.

        The vacuum parameter is used to set the cell size.

        A tuple of exception types can be provided which will be
        caught during a batch of calculations.  Types not specified
        will be considered fatal."""

        dir, path = os.path.split(name)
        self.dir = dir
        self.name = name
        self.vacuum = vacuum
        if exceptions is None:
            exceptions = ()
        self.exceptions = exceptions

    def setup_calculator(self, system, formula):
        """Create a new calculator.

        Default is an EMT calculator.  Most implementations will want to
        override this method."""
        raise NotImplementedError

    def setup_system(self, formula):
        """Create an Atoms object from the given formula.

        By default this will be loaded from the g2 database, setting
        the cell size by means of the molecule test's vacuum parameter."""
        system = molecule(formula)
        system.center(vacuum=self.vacuum)
        return system

    def setup(self, formula):
        """Build calculator and atoms objects.

        This will invoke the setup_calculator and setup_system methods."""
        system = self.setup_system(formula)
        calc = self.setup_calculator(system, formula)
        system.set_calculator(calc)
        return system, calc
        
    def get_filename(self, formula, extension='traj'):
        """Returns the filename for a test result file.

        Default format is <name>.<formula>.traj

        The test may write other files, but this filename is used as a
        flag denoting whether the calculation has been done
        already."""
        return '.'.join([self.name, formula, extension])

    def run(self, formula, system, filename):
        raise NotImplementedError

    def retrieve(self, formula, filename):
        """Retrieve results of previous calculation from file.

        Default implementation returns the total energy.

        This method should be overridden whenever the test method is
        overridden to calculate something else than the total energy."""
        raise NotImplementedError


class EnergyTest:
    def run(self, formula, system, filename):
        """Calculate energy of specified system and save to file."""
        system.get_potential_energy()
         # Won't create .bak file:
        traj = PickleTrajectory(open(filename, 'w'), 'w')
        traj.write(system)
        traj.close()

    def retrieve(self, formula, filename):
        system = read(filename)
        energy = system.get_potential_energy()
        return energy
    
    def calculate_atomization_energies(self, molecular_energies,
                                       atomic_energies):
        atomic_energy_dict = dict(atomic_energies)
        for formula, molecular_energy in molecular_energies:
            try:
                system = molecule(formula)
                atomic = [atomic_energy_dict[s]
                          for s in system.get_chemical_symbols()]            
                atomization_energy = molecular_energy - sum(atomic)
                yield formula, atomization_energy
            except KeyError:
                pass


class BondLengthTest:
    def run(self, formula, system, filename):
        """Calculate bond length of a dimer.

        This will calculate total energies for varying atomic
        separations close to the g2 bond length, allowing
        determination of bond length by fitting.
        """
        if len(system) != 2:
            raise ValueError('Not a dimer')
        traj = PickleTrajectory(open(filename, 'w'), 'w')
        pos = system.positions
        d = np.linalg.norm(pos[1] - pos[0])
        for x in range(-2, 3):
            system.set_distance(0, 1, d * (1.0 + x * 0.02))
            traj.write(system)
        traj.close()
    
    def retrieve(self, formula, filename):
        traj = PickleTrajectory(filename, 'r')
        distances = np.array([np.linalg.norm(a.positions[1] - a.positions[0])
                              for a in traj])
        energies = np.array([a.get_potential_energy() for a in traj])
        polynomial = np.polyfit(distances, energies, 2) # or maybe 3rd order?
        # With 3rd order it is not always obvious which root is right
        pderiv = np.polyder(polynomial, 1)
        d0 = np.roots(pderiv)
        e0 = np.polyval(energies, d0)
        return distances, energies, d0, e0, polynomial
    

class EMTTest(MoleculeTest):
    def setup_calculator(self, system, calculator):
        return EMT()


class EMTEnergyTest(EnergyTest, EMTTest):
    pass


class EMTBondLengthTest(BondLengthTest, EMTTest):
    pass


def main():
    supported_elements = 'Ni, C, Pt, Ag, H, Al, O, N, Au, Pd, Cu'.split(', ')
    formulas = [formula for formula in g1
                if np.all([symbol in supported_elements
                           for symbol
                           in molecule(formula).get_chemical_symbols()])]
    
    atoms = [symbol for symbol in g2_atoms if symbol in supported_elements]
    dimers = [formula for formula in formulas if len(molecule(formula)) == 2]


    name1 = 'testfiles/energy'
    name2 = 'testfiles/bond'
    test1 = BatchTest(EMTEnergyTest(name1, vacuum=3.0))
    test2 = BatchTest(EMTBondLengthTest(name2, vacuum=3.0))

    print 'Energy test'
    print '-----------'
    test1.run(formulas + atoms)

    print
    print 'Bond length test'
    print '----------------'
    test2.run(dimers)

    print
    print 'Atomization energies'
    print '--------------------'
    atomic_energies = dict(test1.collect(atoms))
    molecular_energies = dict(test1.collect(formulas))
    atomization_energies = {}
    for formula, energy in molecular_energies.iteritems():
        system = molecule(formula)
        atomic = [atomic_energies[s] for s in system.get_chemical_symbols()]
        atomization_energy = energy - sum(atomic)
        atomization_energies[formula] = atomization_energy
        print formula.rjust(10), '%.02f' % atomization_energy

    print
    print 'Bond lengths'
    print '------------'
    for formula, (d_i, e_i, d0, e0, poly) in test2.collect(dimers):
        system = molecule(formula)
        bref = np.linalg.norm(system.positions[1] - system.positions[0])
        print formula.rjust(10), '%6.3f' % d0, '  g2ref =', '%2.3f' % bref

        
if __name__ == '__main__':
    main()
