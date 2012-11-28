from ase import Atoms
from ase.io import PickleTrajectory

class Foo(object):
    def __init__(self, value):
        self.value = value
    def __cmp__(self, other):
        return int(self.value - other.value)


if __name__ == '__main__':
    import info  # import ourselves to make info.Foo reachable

    # Create a molecule with an info attribute
    info = dict(creation_date='2011-06-27', 
                chemical_name='Hydrogen',
                # custom classes also works provided that it is
                # imported and pickleable...
                foo=info.Foo(7),  
                )
    molecule = Atoms('H2', positions=[(0., 0., 0.), (0., 0., 1.1)], info=info)
    assert molecule.info == info

    # Copy molecule
    atoms = molecule.copy()
    assert atoms.info == info

    # Save molecule to trajectory
    traj = PickleTrajectory('info.traj', 'w', atoms=molecule)
    traj.write()
    del traj

    # Load molecule from trajectory 
    t = PickleTrajectory('info.traj')
    atoms = t[-1]
    assert atoms.info == info

