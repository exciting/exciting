# -*- coding: utf-8 -*-
"""This module defines I/O routines with CASTEP files.
The key idea is that all function accept or return  atoms objects.
CASTEP specific parameters will be returned through the <atoms>.calc
attribute.
"""

from numpy import sqrt, radians, sin, cos, matrix, array, cross, float32, dot
import ase
from ase.constraints import FixAtoms, FixCartesian
from ase.parallel import paropen
import os


__all__ = [
        'read_castep',
        'read_cell',
        'read_geom',
        'read_param',
        'read_seed',
        'write_cell',
        'write_param',
        ]


def write_cell(filename, atoms, positions_frac=False, castep_cell=None,
    force_write=False):
    """This CASTEP export function write minimal information to
    a .cell file. If the atoms object is a trajectory, it will
    take the last image.
    """
    if atoms is None:
        print("Atoms object not initialized")
        return  False
    if isinstance(atoms, list):
        if len(atoms) > 1:
            atoms = atoms[-1]

    if os.path.isfile(filename) and not force_write:
        print('ase.io.castep.write_param: Set optional argument')
        print('force_write=True to overwrite %s.' % filename)
        return False

    fd = open(filename, 'w')
    fd.write('#######################################################\n')
    fd.write('#CASTEP cell file: %s\n' % filename)
    fd.write('#Created using the Atomic Simulation Environment (ASE)#\n')
    fd.write('#######################################################\n\n')
    fd.write('%BLOCK LATTICE_CART\n')
    cell = matrix(atoms.get_cell())
    for line in atoms.get_cell():
        fd.write('    %.10f %.10f %.10f\n' % tuple(line))
    fd.write('%ENDBLOCK LATTICE_CART\n\n\n')

    if positions_frac:
        keyword = 'POSITIONS_FRAC'
        positions = array(atoms.get_positions() * cell.I)

    else:
        keyword = 'POSITIONS_ABS'
        positions = atoms.get_positions()

    if atoms.get_initial_magnetic_moments().any():
        pos_block = [('%s %8.6f %8.6f %8.6f SPIN=%4.2f' %
            (x, y[0], y[1], y[2], m)) for (x, y, m)
            in zip(atoms.get_chemical_symbols(),
            positions,
            atoms.get_initial_magnetic_moments())]
    else:
        pos_block = [('%s %8.6f %8.6f %8.6f' %
            (x, y[0], y[1], y[2])) for (x, y)
            in zip(atoms.get_chemical_symbols(),
            positions)]

    fd.write('%%BLOCK %s\n' % keyword)
    for line in pos_block:
        fd.write('    %s\n' % line)
    fd.write('%%ENDBLOCK %s\n\n' % keyword)

    # if atoms, has a CASTEP calculator attached, then only
    # write constraints if really necessary
    if hasattr(atoms, 'calc')\
        and hasattr(atoms.calc, 'param')\
        and hasattr(atoms.calc.param, 'task'):
        task = atoms.calc.param.task
        if atoms.calc.param.task.value is None:
            suppress_constraints = True
        elif task.value.lower() not in [
                            'geometryoptimization',
                            'moleculardynamics',
                            'transitionstatesearch',
                            'phonon']:
            suppress_constraints = True
        else:
            suppress_constraints = False
    else:
        suppress_constraints = True

    constraints = atoms.constraints
    if len(constraints) and not suppress_constraints:
        fd.write("%BLOCK IONIC_CONSTRAINTS \n")
        count = 0
        for constr in constraints:
            if not isinstance(constr, FixAtoms)\
                and not isinstance(constr, FixCartesian)\
                and not suppress_constraints:
                print('Warning: you have constraints in your atoms, that are')
                print('         not supported by CASTEP')
                break
            if isinstance(constr, FixAtoms):
                # sorry, for this complicated block
                # reason is that constraint.index can either
                # hold booleans or integers and in both cases
                # it is an numpy array, so no simple comparison works
                for n, val in enumerate(constr.index):
                    if val.dtype.name.startswith('bool'):
                        if not val:
                            continue
                        symbol = atoms.get_chemical_symbols()[n]
                        nis = atoms.calc._get_number_in_species(n)
                    elif val.dtype.name.startswith('int'):
                        symbol = atoms.get_chemical_symbols()[val]
                        nis = atoms.calc._get_number_in_species(val)
                    else:
                        raise UserWarning('Unrecognized index in' + \
                                           ' constraint %s' % constr)
                    fd.write("%6d %3s %3d   1 0 0 \n" % (count + 1,
                                                         symbol,
                                                         nis))
                    fd.write("%6d %3s %3d   0 1 0 \n" % (count + 2,
                                                         symbol,
                                                         nis))
                    fd.write("%6d %3s %3d   0 0 1 \n" % (count + 3,
                                                         symbol,
                                                         nis))
                    count += 3
            elif isinstance(constr, FixCartesian):
                n = constr.a
                symbol = atoms.get_chemical_symbols()[n]
                nis = atoms.calc._get_number_in_species(n)
                fix_cart = - constr.mask + 1
                if fix_cart[0]:
                    count += 1
                    fd.write("%6d %3s %3d   1 0 0 \n" % (count, symbol, nis))
                if fix_cart[1]:
                    count += 1
                    fd.write("%6d %3s %3d   0 1 0 \n" % (count, symbol, nis))
                if fix_cart[2]:
                    count += 1
                    fd.write("%6d %3s %3d   0 0 1 \n" % (count, symbol, nis))
        fd.write("%ENDBLOCK IONIC_CONSTRAINTS \n")

    if castep_cell is None:
        if hasattr(atoms, 'calc') and hasattr(atoms.calc, 'cell'):
            castep_cell = atoms.calc.cell
        else:
            fd.close()
            return True

    for option in castep_cell._options.values():
        if option.value is not None:
            if option.type == 'Block':
                fd.write('%%BLOCK %s\n' % option.keyword.upper())
                fd.write(option.value)
                fd.write('\n%%ENDBLOCK %s\n' % option.keyword.upper())
            else:
                fd.write('%s : %s\n' % (option.keyword.upper(), option.value))

    fd.close()
    return True


def read_cell(filename, _=None):
    """Read a .cell file and return an atoms object.
    Any value found that does not fit the atoms API
    will be stored in the atoms.calc attribute.
    """

    from ase.calculators.castep import Castep
    calc = Castep()

    fileobj = open(filename)
    lines = fileobj.readlines()
    fileobj.close()

    def get_tokens(lines, l):
        """Tokenizes one line of a *cell file."""
        comment_chars = "#!"
        while l < len(lines):
            line = lines[l].strip()
            if len(line) == 0:
                l += 1
                continue
            elif any([line.startswith(comment_char)
                      for comment_char in comment_chars]):
                l += 1
                continue
            else:
                for c in comment_chars:
                    if c in line:
                        icomment = min(line.index(c))
                    else:
                        icomment = len(line) 
                tokens = line[:icomment].split()
                return tokens, l + 1
        tokens = ""
        print("read_cell: Warning - get_tokens has not found any more tokens")
        return tokens, l

    lat = []
    have_lat = False

    pos = []
    spec = []
    constraints = []
    raw_constraints = {}
    have_pos = False
    pos_frac = False

    l = 0
    while l < len(lines):
        tokens, l = get_tokens(lines, l)
        if not tokens:
            continue
        elif tokens[0].upper() == "%BLOCK":
            if tokens[1].upper() == "LATTICE_CART" and not have_lat:
                tokens, l = get_tokens(lines, l)
                if len(tokens) == 1:
                    print('read_cell: Warning - ignoring unit specifier in')
                    print('%BLOCK LATTICE_CART (assuming Angstrom instead)')
                    tokens, l = get_tokens(lines, l)
                for _ in range(3):
                    lat_vec = map(float, tokens[0:3])
                    lat.append(lat_vec)
                    tokens, l = get_tokens(lines, l)
                if tokens[0].upper() != "%ENDBLOCK":
                    print('read_cell: Warning - ignoring more than three')
                    print('lattice vectors in invalid %BLOCK LATTICE_CART')
                    print('%s ...' % tokens[0].upper())
                have_lat = True

            elif tokens[1].upper() == "LATTICE_ABC" and not have_lat:
                tokens, l = get_tokens(lines, l)
                if len(tokens) == 1:
                    print('read_cell: Warning - ignoring unit specifier in')
                    print('%BLOCK LATTICE_ABC (assuming Angstrom instead)')
                    tokens, l = get_tokens(lines, l)
                a, b, c = map(float, tokens[0:3])
                tokens, l = get_tokens(lines, l)
                alpha, beta, gamma = map(lambda phi: radians(float(phi)),
                                                             tokens[0:3])
                tokens, l = get_tokens(lines, l)
                if tokens[0].upper() != "%ENDBLOCK":
                    print('read_cell: Warning - ignoring additional lines in')
                    print('invalid %BLOCK LATTICE_ABC')
                lat_a = [a, 0, 0]
                lat_b = [b * cos(gamma), b * sin(gamma), 0]
                lat_c1 = c * cos(beta)
                lat_c2 = c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma)
                lat_c3 = sqrt(c * c - lat_c1 * lat_c1 - lat_c2 * lat_c2)
                lat_c = [lat_c1, lat_c2, lat_c3]
                lat = [lat_a, lat_b, lat_c]
                have_lat = True

            elif tokens[1].upper() == "POSITIONS_ABS" and not have_pos:
                tokens, l = get_tokens(lines, l)
                if len(tokens) == 1:
                    print('read_cell: Warning - ignoring unit specifier in')
                    print('%BLOCK POSITIONS_ABS(assuming Angstrom instead)')
                    tokens, l = get_tokens(lines, l)
                while len(tokens) == 4:
                    spec.append(tokens[0])
                    pos.append(map(float, tokens[1:4]))
                    tokens, l = get_tokens(lines, l)
                if tokens[0].upper() != "%ENDBLOCK":
                    print('read_cell: Warning - ignoring invalid lines in')
                    print('%%BLOCK POSITIONS_ABS:\n\t %s' % tokens)
                have_pos = True

            elif tokens[1].upper() == "POSITIONS_FRAC" and not have_pos:
                pos_frac = True
                tokens, l = get_tokens(lines, l)
                while len(tokens) == 4:
                    spec.append(tokens[0])
                    pos.append(map(float, tokens[1:4]))
                    tokens, l = get_tokens(lines, l)
                if tokens[0].upper() != "%ENDBLOCK":
                    print('read_cell: Warning - ignoring invalid lines')
                    print('%%BLOCK POSITIONS_FRAC:\n\t %s' % tokens)
                have_pos = True
            elif tokens[1].upper() == 'SPECIES_POT':
                tokens, l = get_tokens(lines, l)
                while tokens and not tokens[0].upper() == '%ENDBLOCK':
                    if len(tokens) == 2:
                        calc.cell.species_pot = tuple(tokens)
                    tokens, l = get_tokens(lines, l)
            elif tokens[1].upper() == 'IONIC_CONSTRAINTS':

                while True:
                    if tokens and tokens[0].upper() == '%ENDBLOCK':
                        break
                    tokens, l = get_tokens(lines, l)
                    if not len(tokens) == 6:
                        continue
                    _, species, nic, x, y, z = tokens
                    nic = int(nic)
                    if (species, nic) not in raw_constraints:
                        raw_constraints[(species, nic)] = []
                    raw_constraints[(species, nic)].append(array(
                                                           [x, y, z]))
            else:
                print('Warning: the keyword %s is not' % tokens[1].upper())
                print('         interpreted in cell files')
                while not tokens[0].upper() == '%ENDBLOCK':
                    tokens, l = get_tokens(lines, l)
                #raise UserWarning
        else:
            key = tokens[0]
            value = ' '.join(tokens[1:])
            try:
                calc.__setattr__(key, value)
            except:
                print("Problem setting calc.cell.%s = %s" % (key, value))
                raise

    if pos_frac:
        atoms = ase.Atoms(
            calculator=calc,
            cell=lat,
            pbc=True,
            scaled_positions=pos,
            symbols=spec,
            )
    else:
        atoms = ase.Atoms(
            calculator=calc,
            cell=lat,
            pbc=True,
            positions=pos,
            symbols=spec,
            )

    fixed_atoms = []
    for (species, nic), value in raw_constraints.iteritems():
        absolute_nr = atoms.calc._get_absolute_number(species, nic)
        if len(value) == 3:
            fixed_atoms.append(absolute_nr)
        elif len(value) == 2:
            constraint = ase.constraints.FixedLine(a=absolute_nr,
                direction=cross(value[0], value[1]))
            constraints.append(constraint)
        elif len(value) == 1:
            constraint = ase.constraints.FixedPlane(a=absolute_nr,
                direction=array(value[0], dtype=float32))
            constraints.append(constraint)
        else:
            print('Error: Found %s statements attached to atoms %s'
    % (len(value), absolute_nr))
    constraints.append(ase.constraints.FixAtoms(fixed_atoms))
    atoms.set_constraint(constraints)

    # needs to go here again to have the constraints in
    # atoms.calc.atoms.constraints as well
    atoms.calc.atoms = atoms
    atoms.calc.push_oldstate()
    return atoms

# this actually does not belong here
# think how one could join this with
# the ase.calculators.castep.Castep.read()
# in the future!


def read_castep(filename, _=-1):
    """Reads a .castep file and returns an atoms  object.
    The calculator information will be stored in the calc attribute.
    If more than one SCF step is found, a list of all steps
    will be stored in the traj attribute.

    Note that the index argument has no effect as of now.
    """
    from ase.calculators.singlepoint import SinglePointCalculator

    fileobj = open(filename)
    lines = fileobj.readlines()
    fileobj.close()
    traj = []
    energy_total = None
    energy_0K = None
    for i, line in enumerate(lines):
        if 'NB est. 0K energy' in line:
            energy_0K = float(line.split()[6])
        elif 'Final energy, E' in line:
            energy_total = float(line.split()[4])
        elif 'Unit Cell' in line:
            cell = [x.split()[0:3] for x in lines[i + 3:i + 6]]
            cell = array([[float(col) for col in row] for row in cell])
        elif 'Cell Contents' in line:
            geom_starts = i
            start_found = False
            for j, jline in enumerate(lines[geom_starts:]):
                if jline.find('xxxxx') > 0 and start_found:
                    geom_stop = j + geom_starts
                    break
                if jline.find('xxxx') > 0 and not start_found:
                    geom_start = j + geom_starts + 4
                    start_found = True
            species = [line.split()[1] for line in lines[geom_start:geom_stop]]
            geom = dot(array([[float(col) for col in line.split()[3:6]]
                for line in lines[geom_start:geom_stop]]), cell)
        elif 'Writing model to' in line:
            atoms = ase.Atoms(
                cell=cell,
                pbc=True,
                positions=geom,
                symbols=''.join(species),
                )
            # take 0K energy where available, else total energy
            if energy_0K:
                energy = energy_0K
            else:
                energy = energy_total
            # generate a minimal single-point calculator
            sp_calc = SinglePointCalculator(atoms=atoms,
                                            energy=energy,
                                            forces=None,
                                            magmoms=None,
                                            stress=None,
                                            )
            atoms.set_calculator(sp_calc)
            traj.append(atoms)

    return traj


def read_param(filename, calc=None):
    """Reads a param file. If an Castep object is passed as the
    second argument, the parameter setings are merged into
    the existing object and returned. Otherwise a new Castep()
    calculator instance gets created and returned.
    Parameters:
        filename: the .param file. Only opens reading
        calc: [Optional] calculator object to hang
            parameters onto
    """
    if calc is None:
        calc = ase.calculators.castep.Castep(check_castep_version=False)
    calc.merge_param(filename)
    return calc


def write_param(filename, param, check_checkfile=False,
                                 force_write=False,
                                 interface_options=None):
    """Writes a CastepParam object to a CASTEP .param file

    Parameters:
        filename: the location of the file to write to. If it
        exists it will be overwritten without warning. If it
        doesn't it will be created.
        param: a CastepParam instance
        check_checkfile : if set to True, write_param will
            only write continuation or reuse statement
            if a restart file exists in the same directory
    """
    if os.path.isfile(filename) and not force_write:
        print('ase.io.castep.write_param: Set optional argument')
        print('force_write=True to overwrite %s.' % filename)
        return False

    out = paropen(filename, 'w')
    out.write('#######################################################\n')
    out.write('#CASTEP param file: %s\n' % filename)
    out.write('#Created using the Atomic Simulation Environment (ASE)#\n')
    if interface_options is not None:
        out.write('# Internal settings of the calculator\n')
        out.write('# This can be switched off by settings\n')
        out.write('# calc._export_settings = False\n')
        out.write('# If stated, this will be automatically processed\n')
        out.write('# by ase.io.castep.read_seed()\n')
        for option, value in sorted(interface_options.iteritems()):
            out.write('# ASE_INTERFACE %s : %s\n' % (option, value))
    out.write('#######################################################\n\n')
    for keyword, opt in sorted(param._options.iteritems()):
        if opt.type == 'Defined':
            if opt.value is not None:
                out.write('%s\n' % (option))
        elif opt.value is not None:
            if keyword in ['continuation', 'reuse'] and check_checkfile:
                if opt.value == 'default':
                    if not os.path.exists('%s.%s'\
                        % (os.path.splitext(filename)[0], 'check')):
                        continue
                elif not os.path.exists(opt.value):
                    continue
            out.write('%s : %s\n'
                % (keyword, opt.value))
    out.close()


def read_geom(filename, _=-1):
    """Reads a .geom file produced by the CASTEP GeometryOptimization task and
    returns an atoms  object.
    The information about total free energy and forces of each atom for every
    relaxation step will be stored for further analysis especially in a
    single-point calculator.
    Note that everything in the .geom file is in atomic units, which has
    been conversed to commonly used unit angstrom(length) and eV (energy).

    Note that the index argument has no effect as of now.

    Contribution by Wei-Bing Zhang. Thanks!
    """
    from ase.calculators.singlepoint import SinglePointCalculator
    fileobj = open(filename)
    txt = fileobj.readlines()
    fileobj.close()
    traj = []

    # Source: CODATA2002, used by default
    # in CASTEP 5.01
    # but check with your version in case of error
    # ase.units is based on CODATA1986/
    # here we hard-code from http://physics.nist.gov/cuu/Document/all_2002.pdf
    Hartree = 27.211384565719481
    Bohr = 0.5291772108

    print('N.B.: Energy in .geom file is not 0K extrapolated.')
    for i, line in enumerate(txt):
        if line.find("<-- E") > 0:
            start_found = True
            energy = float(line.split()[0]) * Hartree
            cell = [x.split()[0:3] for x in txt[i + 1:i + 4]]
            cell = array([[float(col) * Bohr for col in row] for row in
                cell])
        if line.find('<-- R') > 0 and start_found:
            start_found = False
            geom_start = i
            for i, line in enumerate(txt[geom_start:]):
                if line.find('<-- F') > 0:
                    geom_stop = i + geom_start
                    break
            species = [line.split()[0] for line in
                txt[geom_start:geom_stop]]
            geom = array([[float(col) * Bohr for col in
                line.split()[2:5]] for line in txt[geom_start:geom_stop]])
            forces = array([[float(col) * Hartree / Bohr for col in
                line.split()[2:5]] for line in
                    txt[geom_stop:geom_stop + (geom_stop - geom_start)]])
            image = ase.Atoms(species, geom, cell=cell, pbc=True)
            image.set_calculator(SinglePointCalculator(energy, forces, None,
                None, image))
            traj.append(image)

    return traj


def read_seed(seed, new_seed=None, ignore_internal_keys=False):
    """A wrapper around the CASTEP Calculator in conjunction with
    read_cell and read_param. Basically this can be used to reuse
    a previous calculation which results in a triple of
    cell/param/castep file. The label of the calculation if pre-
    fixed with cop_of_ and everything else will be recycled as
    much as possible from the addressed calculation.
    """

    directory = os.path.abspath(os.path.dirname(seed))
    seed = os.path.basename(seed)

    paramfile = os.path.join(directory, '%s.param' % seed)
    cellfile = os.path.join(directory, '%s.cell' % seed)
    castepfile = os.path.join(directory, '%s.castep' % seed)

    atoms = read_cell(cellfile)
    atoms.calc._directory = directory
    atoms.calc._rename_existing_dir = False
    atoms.calc._castep_pp_path = directory
    atoms.calc.merge_param(paramfile,
                           ignore_internal_keys=ignore_internal_keys)
    if new_seed is None:
        atoms.calc._label = 'copy_of_%s' % seed
    else:
        atoms.calc._label = str(new_seed)
    if os.path.isfile(castepfile):
        # _set_atoms needs to be True here
        # but we set it right back to False
        atoms.calc._set_atoms = True
        atoms.calc.read(castepfile)
        atoms.calc._set_atoms = False
        # sync the top-level object with the
        # one attached to the calculator
        atoms = atoms.calc.atoms
    else:
        print('Corresponding CASTEP not found.')
    atoms.calc.push_oldstate()

    return atoms
