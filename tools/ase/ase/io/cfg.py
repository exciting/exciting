import numpy as np

import ase

from ase.parallel import paropen


cfg_default_fields = np.array( [ 'positions', 'momenta', 'numbers', 'magmoms' ] )

def write_cfg(f, a):
    """Write atomic configuration to a CFG-file (native AtomEye format).
       See: http://mt.seas.upenn.edu/Archive/Graphics/A/
    """
    if isinstance(f, str):
        f = paropen(f, 'w')

    f.write('Number of particles = %i\n' % len(a))
    f.write('A = 1.0 Angstrom\n')
    cell = a.get_cell()
    for i in range(3):
        for j in range(3):
            f.write('H0(%1.1i,%1.1i) = %f A\n' % ( i + 1, j + 1, cell[i, j] ))

    entry_count  = 3
    for x in a.arrays.keys():
        if not x in cfg_default_fields:
            if len(a.get_array(x).shape) == 1:
                entry_count += 1
            else:
                entry_count += a.get_array(x).shape[1]

    vels = a.get_velocities()
    if type(vels) == np.ndarray:
        entry_count += 3
    else:
        f.write('.NO_VELOCITY.\n')

    f.write('entry_count = %i\n' % entry_count)

    i = 0
    for name, aux in a.arrays.iteritems():
        if not name in cfg_default_fields:
            if len(aux.shape) == 1:
                f.write('auxiliary[%i] = %s [a.u.]\n' % ( i, name ))
                i += 1
            else:
                for j in range(aux.shape[1]):
                    f.write('auxiliary[%i] = %s_%1.1i [a.u.]\n' % ( i, name, j ))
                    i += 1

    # Distinct elements
    spos = a.get_scaled_positions()
    for i in a:
        el  = i.symbol

        f.write('%f\n' % ase.data.atomic_masses[ase.data.chemical_symbols.index(el)])
        f.write('%s\n' % el)

        x, y, z = spos[i.index, :]
        s =  '%e %e %e ' % ( x, y, z )

        if type(vels) == np.ndarray:
            vx, vy, vz  = vels[i.index, :]
            s  = s + ' %e %e %e ' % ( vx, vy, vz )

        for name, aux in a.arrays.iteritems():
            if not name in cfg_default_fields:
                if len(aux.shape) == 1:
                    s += ' %e' % aux[i.index]
                else:
                    s += ( aux.shape[1]*' %e' ) % tuple(aux[i.index].tolist())
                         
        f.write('%s\n' % s)


default_color = {
    'H': [ 0.800, 0.800, 0.800 ],
    'C': [ 0.350, 0.350, 0.350 ],
    'O': [ 0.800, 0.200, 0.200 ]
    }

default_radius = {
    'H': 0.435,
    'C': 0.655,
    'O': 0.730
    }
       


def write_clr(f, atoms):
    """Write extra color and radius code to a CLR-file (for use with AtomEye).
       Hit F12 in AtomEye to use.
       See: http://mt.seas.upenn.edu/Archive/Graphics/A/
    """
    color   = None
    radius  = None
    if atoms.has('color'):
        color  = atoms.get_array('color')
    if atoms.has('radius'):
        radius  = atoms.get_array('radius')

    if color is None:
        color  = np.zeros([len(atoms),3], dtype=float)
        for a in atoms:
            color[a.index, :]  = default_color[a.symbol]

    if radius is None:
        radius  = np.zeros(len(atoms), dtype=float)
        for a in atoms:
            radius[a.index]  = default_radius[a.symbol]

    radius.shape = (-1, 1)

    if isinstance(f, str):
        f = paropen(f, 'w')   
    for c1, c2, c3, r in np.append(color, radius, axis=1):
        f.write('%f %f %f %f\n' % ( c1, c2, c3, r ))


###

def read_key_val(f):
    if isinstance(f, str):
        l = f
    else:
        l = f.readline()
    s = l.split('=')
    if len(s) != 2:
        raise RuntimeError("Line '%s' is not of the form 'key = value'." % l[:-1])
    return ( s[0].strip(), s[1].strip() )


def read_str_key(f, key, key2=None):
    in_key, val  = read_key_val(f)
    if key2 is None:
        if key.upper() != in_key.upper():
            raise RuntimeError("Key '%s' expected, '%s' found." % ( key, in_key ))
    else:
        if key.upper() != in_key.upper() and key2.upper() != in_key.upper():
            raise RuntimeError("Key '%s' or '%s' expected, '%s' found." % ( key, key2, in_key ))
    return val


def read_int_key(f, key):
    vals  = read_str_key(f, key).split()
    # Ignore units
    return int(vals[0])


def read_float_key(f, key):
    vals  = read_str_key(f, key).split()
    # Ignore units
    return float(vals[0])


###

def read_cfg(f):
    """Read atomic configuration from a CFG-file (native AtomEye format).
       See: http://mt.seas.upenn.edu/Archive/Graphics/A/
    """
    if isinstance(f, str):
        f  = open(f)

    nat   = read_int_key(f, 'Number of particles')
    unit  = read_float_key(f, 'A')

    cell  = np.zeros( [ 3, 3 ] )
    for i in range(3):
        for j in range(3):
            cell[i, j]  = read_float_key(f, 'H0(%1.1i,%1.1i)' % (i + 1, j + 1))

    l     = f.readline()
    vels  = None
    if l.strip() == '.NO_VELOCITY.':
        l     = f.readline()
    else:
        vels  = np.zeros( [ nat, 3 ] )

    naux      = read_int_key(l, 'entry_count') - 3
    if vels is not None:
        naux -= 3
    aux       = np.zeros( [ nat, naux ] )

    auxstrs   = [ ]
    for i in range(naux):
        s  = read_str_key(f, 'auxiliary[%1.1i]' % i, 'auxiliary[%2.2i]' % i)
        auxstrs += [ s[:s.find('[')].strip() ]

    spos     = np.zeros( [ nat, 3 ] )
    masses   = np.zeros( nat )
    syms     = [ '' for i in range(nat) ]

    i  = 0
    s  = f.readline().split()
    while l:
        mass   = float(s[0])
        sym    = f.readline().strip()

        l      = f.readline()
        s      = l.split()

        while l and len(s) > 1:
            masses[i]  = mass
            syms[i]    = sym
            props      = [ float(x) for x in s ]

            spos[i, :] = props[0:3]
            off        = 3
            if vels is not None:
                off         = 6
                vels[i, :]  = props[3:6]

            aux[i, :]  = props[off:]

            i  += 1

            l  = f.readline()
            if l:
                s  = l.split()
                
    if vels is None:
        a = ase.Atoms(
            symbols           = syms,
            masses            = masses,
            scaled_positions  = spos,
            cell              = cell,
            pbc               = True
        )
    else:
        a = ase.Atoms(
            symbols           = syms,
            masses            = masses,
            scaled_positions  = spos,
            momenta           = masses.reshape(-1,1)*vels,
            cell              = cell,
            pbc               = True
        )

    i  = 0
    while i < naux:
        auxstr  = auxstrs[i]

        if auxstr[-2:] == '_x':
            a.set_array(auxstr[:-2], aux[:, i:i+3])

            i  += 3
        else:
            a.set_array(auxstr, aux[:, i])

            i  += 1

    return a
