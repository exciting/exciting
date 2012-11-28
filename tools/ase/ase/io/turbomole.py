from ase.atoms import Atoms
from ase.units import Bohr


def read_turbomole(filename='coord'):
    """Method to read turbomole coord file
    
    coords in bohr, atom types in lowercase, format:
    $coord
    x y z atomtype 
    x y z atomtype f
    $end
    Above 'f' means a fixed atom.
    """
    from ase import Atoms, Atom
    from ase.constraints import FixAtoms

    if isinstance(filename, str):
        f = open(filename)

    lines = f.readlines()
    atoms_pos = []
    atom_symbols = []
    myconstraints=[]
    
    # find $coord section;
    # does not necessarily have to be the first $<something> in file...
    for i, l in enumerate(lines):
        if l.strip().startswith('$coord'):
            start = i
            break
    for line in lines[start+1:]:
        if line.startswith('$'): # start of new section
            break
        else:
            x, y, z, symbolraw = line.split()[:4]
            symbolshort=symbolraw.strip()
            symbol=symbolshort[0].upper()+symbolshort[1:].lower()
            #print symbol
            atom_symbols.append(symbol)
            atoms_pos.append([float(x)*Bohr, float(y)*Bohr, float(z)*Bohr])
            cols = line.split()
            if (len(cols) == 5):
                fixedstr = line.split()[4].strip()
                if (fixedstr == "f"):
                    myconstraints.append(True)
                else:
                    myconstraints.append(False)
            else:
                myconstraints.append(False)
            
    if type(filename) == str:
        f.close()

    atoms = Atoms(positions = atoms_pos, symbols = atom_symbols, pbc = False)
    c = FixAtoms(mask = myconstraints)
    atoms.set_constraint(c)
    #print c
    

    return atoms

def read_turbomole_gradient(filename='gradient', index=-1):
    """ Method to read turbomole gradient file """

    if isinstance(filename, str):
        f = open(filename)

    # read entire file
    lines = [x.strip() for x in f.readlines()]

    # find $grad section
    start = end = -1
    for i, line in enumerate(lines):
        if not line.startswith('$'):
            continue
        if line.split()[0] == '$grad':
            start = i
        elif start >= 0:
            end = i
            break

    if end <= start:
        raise RuntimeError('File %s does not contain a valid \'$grad\' section' % (filename))

    def formatError():
        raise RuntimeError('Data format in file %s does not correspond to known Turbomole gradient format' % (filename))


    # trim lines to $grad
    del lines[:start+1]
    del lines[end-1-start:]

    # Interpret $grad section
    from ase import Atoms, Atom
    from ase.calculators.singlepoint import SinglePointCalculator
    from ase.units import Bohr
    images = []
    while len(lines): # loop over optimization cycles
        # header line
        # cycle =      1    SCF energy =     -267.6666811409   |dE/dxyz| =  0.157112
        fields = lines[0].split('=')
        try:
            cycle = int(fields[1].split()[0])
            energy = float(fields[2].split()[0])
            gradient = float(fields[3].split()[0])
        except (IndexError, ValueError):
            formatError()
        
        # coordinates/gradient
        atoms = Atoms()
        forces = []
        for line in lines[1:]:
            fields = line.split()
            if len(fields) == 4: # coordinates
                # 0.00000000000000      0.00000000000000      0.00000000000000      c
                try:
                    symbol = fields[3].lower().capitalize()
                    position = tuple([bohr2angstrom(float(x)) for x in fields[0:3] ])
                except ValueError:
                    formatError()
                atoms.append(Atom(symbol, position))
            elif len(fields) == 3: # gradients
                #  -.51654903354681D-07  -.51654903206651D-07  0.51654903169644D-07
                try:
                    grad = [float(x.replace('D', 'E')) * Bohr for x in fields[0:3] ]
                except ValueError:
                    formatError()
                forces.append(grad)
            else: # next cycle
                break

        # calculator
        calc = SinglePointCalculator(energy, forces, None, None, atoms)
        atoms.set_calculator(calc)

        # save frame
        images.append(atoms)

        # delete this frame from data to be handled
        del lines[:2*len(atoms)+1]

    return images[index]


def write_turbomole(filename, atoms):
    """Method to write turbomole coord file
    """

    import numpy as np
    from ase.constraints import FixAtoms

    if isinstance(filename, str):
        f = open(filename, 'w')
    else: # Assume it's a 'file-like object'
        f = filename

    coord = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    printfixed = False

    if atoms.constraints:
        for constr in atoms.constraints:
            if isinstance(constr, FixAtoms):
                fix_index=constr.index
                printfixed=True
    #print sflags
        
    if (printfixed):
        fix_str=[]
        for i in fix_index:
            if i == 1:
                fix_str.append("f")
            else:
                fix_str.append(" ")


    f.write("$coord\n")
    if (printfixed):
        for (x, y, z), s, fix in zip(coord,symbols,fix_str):
            f.write('%20.14f  %20.14f  %20.14f      %2s  %2s \n' 
                    % (x/Bohr, y/Bohr, z/Bohr, s.lower(), fix))

    else:
        for (x, y, z), s in zip(coord,symbols):
            f.write('%20.14f  %20.14f  %20.14f      %2s \n' 
                    % (x/Bohr, y/Bohr, z/Bohr, s.lower()))
    f.write("$end\n")
    
