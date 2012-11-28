
def read_aims(filename):
    """Import FHI-aims geometry type files.

    Reads unitcell, atom positions and constraints from
    a geometry.in file.
    """

    from ase import Atoms
    from ase.constraints import FixAtoms, FixCartesian
    import numpy as np

    atoms = Atoms()
    fd = open(filename, 'r')
    lines = fd.readlines()
    fd.close()
    positions = []
    cell = []
    symbols = []
    fix = []
    fix_cart = []
    xyz = np.array([0, 0, 0])
    i = -1
    n_periodic = -1
    periodic = np.array([False, False, False])
    for n, line in enumerate(lines):
        inp = line.split()
        if inp == []:
            continue
        if inp[0] == 'atom':
            if xyz.all():
                fix.append(i)
            elif xyz.any():
                fix_cart.append(FixCartesian(i, xyz))
            floatvect = float(inp[1]), float(inp[2]), float(inp[3])
            positions.append(floatvect)
            symbols.append(inp[-1])
            i += 1
            xyz = np.array([0, 0, 0])
        elif inp[0] == 'lattice_vector':
            floatvect = float(inp[1]), float(inp[2]), float(inp[3])
            cell.append(floatvect)
            n_periodic = n_periodic + 1
            periodic[n_periodic] = True
        if inp[0] == 'constrain_relaxation':
            if inp[1] == '.true.':
                fix.append(i)
            elif inp[1] == 'x':
                xyz[0] = 1
            elif inp[1] == 'y':
                xyz[1] = 1
            elif inp[1] == 'z':
                xyz[2] = 1
    if xyz.all():
        fix.append(i)
    elif xyz.any():
        fix_cart.append(FixCartesian(i, xyz))
    atoms = Atoms(symbols, positions)
    if periodic.all():
        atoms.set_cell(cell)
        atoms.set_pbc(periodic)
    if len(fix):
        atoms.set_constraint([FixAtoms(indices=fix)]+fix_cart)
    else:
        atoms.set_constraint(fix_cart)
    return atoms

def write_aims(filename, atoms):
    """Method to write FHI-aims geometry files.

    Writes the atoms positions and constraints (only FixAtoms is
    supported at the moment). 
    """

    from ase.constraints import FixAtoms, FixCartesian
    import numpy as np

    if isinstance(atoms, (list, tuple)):
        if len(atoms) > 1:
            raise RuntimeError("Don't know how to save more than "+
                               "one image to FHI-aims input")
        else:
            atoms = atoms[0]

    fd = open(filename, 'w')
    fd.write('#=======================================================\n')
    fd.write('#FHI-aims file: '+filename+'\n')
    fd.write('#Created using the Atomic Simulation Environment (ASE)\n')
    fd.write('#=======================================================\n')
    i = 0
    if atoms.get_pbc().any():
        for n, vector in enumerate(atoms.get_cell()):
            fd.write('lattice_vector ')
            for i in range(3):
                fd.write('%16.16f ' % vector[i])
            fd.write('\n')
    fix_cart = np.zeros([len(atoms),3]) 

    if atoms.constraints:
        for constr in atoms.constraints:
            if isinstance(constr, FixAtoms):
                fix_cart[constr.index] = [1,1,1]
            elif isinstance(constr, FixCartesian):
                fix_cart[constr.a] = -constr.mask+1

    for i, atom in enumerate(atoms):
        fd.write('atom ')
        for pos in atom.position:
            fd.write('%16.16f ' % pos)
        fd.write(atom.symbol)
        fd.write('\n')
# (1) all coords are constrained:
        if fix_cart[i].all():
            fd.write('constrain_relaxation .true.\n')
# (2) some coords are constrained:
        elif fix_cart[i].any():
            xyz = fix_cart[i]
            for n in range(3):
                if xyz[n]:
                    fd.write('constrain_relaxation %s\n' % 'xyz'[n])
        if atom.charge:
            fd.write('initial_charge %16.6f\n' % atom.charge)
        if atom.magmom:
            fd.write('initial_moment %16.6f\n' % atom.magmom)
# except KeyError:
#     continue

def read_energy(filename):
    for line in open(filename, 'r'):
        if line.startswith('  | Total energy corrected'):
            E = float(line.split()[-2])
    return E

def read_aims_output(filename, index = -1):
    """  Import FHI-aims output files with all data available, i.e. relaxations, 
    MD information, force information etc etc etc. """
    from ase import Atoms, Atom 
    from ase.calculators.singlepoint import SinglePointCalculator
    from ase.units import Ang, fs
    from ase.constraints import FixAtoms, FixCartesian
    molecular_dynamics = False
    fd = open(filename, 'r')
    cell = []
    images = []
    fix = []
    fix_cart = []    
    n_periodic = -1
    f = None
    pbc = False
    found_aims_calculator = False
    v_unit = Ang/(1000.0*fs)
    while True:
        line = fd.readline()
        if not line:
            break
        if "List of parameters used to initialize the calculator:" in line:
            fd.readline()
            calc = read_aims_calculator(fd)
            calc.out = filename
            found_aims_calculator = True
        if "Number of atoms" in line:
            inp = line.split()
            n_atoms = int(inp[5])
        if "| Unit cell:" in line:
            if not pbc:
                pbc = True
                for i in range(3):
                    inp = fd.readline().split()
                    cell.append([inp[1],inp[2],inp[3]])
        if "Found relaxation constraint for atom" in line:
            xyz = [0, 0, 0]
            ind = int(line.split()[5][:-1])-1
            if "All coordinates fixed" in line:
                if ind not in fix:
                    fix.append(ind)
            if "coordinate fixed" in line:
                coord = line.split()[6]
                constr_ind = 0
                if coord == 'x':
                    xyz[0] = 1
                elif coord == 'y':
                    xyz[1] = 1
                elif coord == 'z':
                    xyz[2] = 1
                keep = True
                for n,c in enumerate(fix_cart):
                    if ind == c.a:
                        keep = False
                        constr_ind = n
                if keep:
                    fix_cart.append(FixCartesian(ind, xyz))
                else:
                    fix_cart[n].mask[xyz.index(1)] = 0
        if "Atomic structure:" in line and not molecular_dynamics:
            fd.readline()
            atoms = Atoms()
            for i in range(n_atoms):
                inp = fd.readline().split()
                atoms.append(Atom(inp[3],(inp[4],inp[5],inp[6])))
        if "Complete information for previous time-step:" in line:
            molecular_dynamics = True
        if "Updated atomic structure:" in line and not molecular_dynamics:
            fd.readline()
            atoms = Atoms()
            velocities = []
            for i in range(n_atoms):
                inp = fd.readline().split()
                if 'lattice_vector' in inp[0]:
                    cell = []
                    for i in range(3):
                        cell += [[float(inp[1]),float(inp[2]),float(inp[3])]]
                        inp = fd.readline().split()
                    atoms.set_cell(cell)
                    inp = fd.readline().split()
                atoms.append(Atom(inp[4],(inp[1],inp[2],inp[3])))  
                if molecular_dynamics:
                    inp = fd.readline().split()
        if "Atomic structure (and velocities)" in line:
            fd.readline()
            atoms = Atoms()
            velocities = []
            for i in range(n_atoms):
                inp = fd.readline().split()
                atoms.append(Atom(inp[4],(inp[1],inp[2],inp[3])))  
                inp = fd.readline().split()
                velocities += [[float(inp[1])*v_unit,float(inp[2])*v_unit,float(inp[3])*v_unit]]
            atoms.set_velocities(velocities)
            if len(fix):
                atoms.set_constraint([FixAtoms(indices=fix)]+fix_cart)
            else:
                atoms.set_constraint(fix_cart)
            images.append(atoms)
        if "Total atomic forces" in line:
            f = []
            for i in range(n_atoms):
                inp = fd.readline().split()
                f.append([float(inp[2]),float(inp[3]),float(inp[4])])
            if not found_aims_calculator:
                e = images[-1].get_potential_energy()
                images[-1].set_calculator(SinglePointCalculator(e,f,None,None,atoms))
            e = None
            f = None
        if "Total energy corrected" in line:
            e = float(line.split()[5])
            if pbc:
                atoms.set_cell(cell)
                atoms.pbc = True
            if not found_aims_calculator:
                atoms.set_calculator(SinglePointCalculator(e,None,None,None,atoms))
            if not molecular_dynamics: 
                if len(fix):
                    atoms.set_constraint([FixAtoms(indices=fix)]+fix_cart)
                else:
                    atoms.set_constraint(fix_cart)
                images.append(atoms)
            e = None
            if found_aims_calculator:
                calc.set_results(images[-1])
                images[-1].set_calculator(calc)
    fd.close()
    if molecular_dynamics:
        images = images[1:]

    # return requested images, code borrowed from ase/io/trajectory.py
    if isinstance(index, int):
        return images[index]
    else:
        step = index.step or 1
        if step > 0:
            start = index.start or 0
            if start < 0:
                start += len(images)
            stop = index.stop or len(images)
            if stop < 0:
                stop += len(images)
        else:
            if index.start is None:
                start = len(images) - 1
            else:
                start = index.start
                if start < 0:
                    start += len(images)
            if index.stop is None:
                stop = -1
            else:
                stop = index.stop
                if stop < 0:
                    stop += len(images)
        return [images[i] for i in range(start, stop, step)]

def read_aims_calculator(file):
    """  found instructions for building an FHI-aims calculator in the output file, 
    read its specifications and return it. """
    from ase.calculators.aims import Aims
    calc = Aims()
    while True:
        line = file.readline()
        if "=======================================================" in line:
            break
        else:
            args = line.split()
            key = '#'
            if len(args) > 0:
                key = args[0]
            if key == '#':
                comment = True   
            elif calc.float_params.has_key(key):
                calc.float_params[key] = float(args[1])
            elif calc.exp_params.has_key(key):
                calc.exp_params[key] = float(args[1])
            elif calc.string_params.has_key(key):
                calc.string_params[key] = args[1]
                if len(args) > 2:
                    for s in args[2:]:
                        calc.string_params[key] += " "+s
            elif calc.int_params.has_key(key):
                calc.int_params[key] = int(args[1])
            elif calc.bool_params.has_key(key):
                try:
                    calc.bool_params[key] = bool(args[1])
                except:
                    if key == 'vdw_correction_hirshfeld':
                        calc.bool_params[key] = True
            elif calc.list_params.has_key(key):
                if key == 'output':
                    # build output string from args:
                    out_option = ''
                    for arg in args[1:]:
                        out_option +=str(arg)+' '
                    if calc.list_params['output'] is not None:
                        calc.list_params['output'] += [out_option]
                    else:
                        calc.list_params['output'] = [out_option]
                else:
                    calc.list_params[key] = list(args[1:])
            elif '#' in key or calc.input_parameters.has_key(key):
                key = key[1:]
                if calc.input_parameters.has_key(key):
                    calc.input_parameters[key] = args[1]
                    if len(args) > 2: 
                        for s in args[2:]:
                            calc.input_parameters[key] += " "+s                
            else:
                raise TypeError('FHI-aims keyword not defined in ASE: ' + key + '. Please check.')
    return calc
