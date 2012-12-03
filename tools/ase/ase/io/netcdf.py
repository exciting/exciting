"""Read and write ASE2's netCDF trajectory files."""

from ase.io.pupynere import NetCDFFile
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator


def read_netcdf(filename, index=-1):
    nc = NetCDFFile(filename)
    dims = nc.dimensions
    vars = nc.variables

    positions = vars['CartesianPositions']
    numbers = vars['AtomicNumbers'][:]
    pbc = vars['BoundaryConditions'][:]
    cell = vars['UnitCell']
    tags = vars['Tags'][:]
    if not tags.any():
        tags = None
    magmoms = vars['MagneticMoments'][:]
    if not magmoms.any():
        magmoms = None

    nimages = positions.shape[0]

    attach_calculator = False
    if 'PotentialEnergy' in vars:
        energy = vars['PotentialEnergy']
        attach_calculator = True
    else:
        energy = nimages * [None]

    if 'CartesianForces' in vars:
        forces = vars['CartesianForces']
        attach_calculator = True
    else:
        forces = nimages * [None]

    if 'Stress' in vars:
        stress = vars['Stress']
        attach_calculator = True
    else:
        stress = nimages * [None]

    if isinstance(index, int):
        indices = [index]
    else:
        indices = range(nimages)[index]

    images = []
    for i in indices:
        atoms = Atoms(positions=positions[i],
                      numbers=numbers,
                      cell=cell[i],
                      pbc=pbc,
                      tags=tags, magmoms=magmoms)

        if attach_calculator:
            calc = SinglePointCalculator(energy[i], forces[i], stress[i],
                                         None, atoms) ### Fixme magmoms
            atoms.set_calculator(calc)
            
        images.append(atoms)
        
    if isinstance(index, int):
        return images[0]
    else:
        return images


class LOA:
    def __init__(self, images):
        self.set_atoms(images[0])

    def __len__(self):
        return len(self.atoms)
    
    def set_atoms(self, atoms):
        self.atoms = atoms
        
    def GetPotentialEnergy(self):
        return self.atoms.get_potential_energy()

    def GetCartesianForces(self):
        return self.atoms.get_forces()

    def GetUnitCell(self):
        return self.atoms.get_cell()

    def GetAtomicNumbers(self):
        return self.atoms.get_atomic_numbers()

    def GetCartesianPositions(self):
        return self.atoms.get_positions()

    def GetBoundaryConditions(self):
        return self.atoms.get_pbc()
    

def write_netcdf(filename, images):
    from ASE.Trajectories.NetCDFTrajectory import NetCDFTrajectory

    if not isinstance(images, (list, tuple)):
        images = [images]
        
    loa = LOA(images)
    traj = NetCDFTrajectory(filename, loa)
    for atoms in images:
        loa.set_atoms(atoms)
        traj.Update()
    traj.Close()

