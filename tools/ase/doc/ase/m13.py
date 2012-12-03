from math import sqrt
from ase.cluster.cubic import FaceCenteredCubic
from ase.tasks.molecule import MoleculeTask
from ase.data import covalent_radii, atomic_numbers

class M13Task(MoleculeTask):
    taskname = 'm13'
    def build_system(self, name):
        if self.bond_length is None:
            b = 2 * covalent_radii[atomic_numbers[name]]
        else:
             b = self.bond_length

        return FaceCenteredCubic(name, [(1, 0, 0)], [1],
                                 latticeconstant=b * sqrt(2))

task = M13Task()

