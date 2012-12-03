
class Calculator:
    def __init__(self):
        return

    def set_atoms(self, atoms):
        self.atoms = atoms.copy()

    def get_atoms(self):
        atoms = self.atoms.copy()
        atoms.set_calculator(self)
        return atoms

    def get_name(self):
        """Return the name of the calculator (string).  """
        raise NotImplementedError

    def get_version(self):
        """Return the version of the calculator (string).  """
        raise NotImplementedError

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)
        if force_consistent:
            return self.energy_free
        else:
            return self.energy_zero

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces

    def get_stress(self, atoms):
        self.update(atoms)
        return self.stress

    def initialize(self, atoms):
        """Prepare the input files required to
        start the program (calculator).  """
        raise NotImplementedError

    def read(self, atoms):
        self.positions = atoms.get_positions()
        self.energy_free, self.energy_zero = self.read_energy()
        self.forces = self.read_forces(atoms)
        self.dipole = self.read_dipole()
        self.fermi = self.read_fermi()
        self.atoms = atoms.copy()
        try:
            self.nbands = self.read_nbands()
        except NotImplementedError:
            do_nothing = True
        except AttributeError:
            do_nothing = True
        try:
            self.stress = self.read_stress()
        except NotImplementedError:
            do_nothing = True
        return
