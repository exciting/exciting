import numpy as np

from ase.calculators.singlepoint import SinglePointCalculator
from ase.atom import Atom
from ase.atoms import Atoms


def read_dacapo_text(fileobj):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    lines = fileobj.readlines()
    i = lines.index(' Structure:             A1           A2            A3\n')
    cell = np.array([[float(w) for w in line.split()[2:5]]
                     for line in lines[i + 1:i + 4]]).transpose()
    i = lines.index(' Structure:  >>         Ionic positions/velocities ' +
                    'in cartesian coordinates       <<\n')
    atoms = []
    for line in lines[i + 4:]:
        words = line.split()
        if len(words) != 9:
            break
        Z, x, y, z = words[2:6]
        atoms.append(Atom(int(Z), [float(x), float(y), float(z)]))

    atoms = Atoms(atoms, cell=cell.tolist())

    try:
        i = lines.index(
            ' DFT:  CPU time                           Total energy\n')
    except ValueError:
        pass
    else:
        column = lines[i + 3].split().index('selfcons') - 1
        try:
            i2 = lines.index(' ANALYSIS PART OF CODE\n', i)
        except ValueError:
            pass
        else:
            while i2 > i:
                if lines[i2].startswith(' DFT:'):
                    break
                i2 -= 1
            energy = float(lines[i2].split()[column])
            atoms.set_calculator(SinglePointCalculator(energy, None, None,
                                                       None, atoms))

    return atoms



def read_dacapo(filename):
    from ase.io.pupynere import NetCDFFile

    nc = NetCDFFile(filename)
    dims = nc.dimensions
    vars = nc.variables

    cell = vars['UnitCell'][-1]
    try:
        magmoms = vars['InitialAtomicMagneticMoment'][:]
    except KeyError:
        magmoms = None
    try:
        tags = vars['AtomTags'][:]
    except KeyError:
        tags = None
    atoms = Atoms(scaled_positions=vars['DynamicAtomPositions'][-1],
                  symbols=[(a + b).strip() 
                           for a, b in vars['DynamicAtomSpecies'][:]],
                  cell=cell,
                  magmoms=magmoms,
                  tags=tags,
                  pbc=True)

    try:
        energy = vars['TotalEnergy'][-1] 
        force = vars['DynamicAtomForces'][-1] 
    except KeyError: 
        energy = None 
        force = None 
    calc = SinglePointCalculator(energy,force,None, None, atoms)  ### Fixme magmoms
    atoms.set_calculator(calc)
        
    return atoms
