def get_ionization_energy(name, vertical=True):
    "Deprecated, use ase.data.g2.get_ionization_energy instead."
    import warnings
    warnings.warn('ase.data.molecules.get_ionization_energy is deprecated. '
                  ' Please use from ase.data.g2 import get_ionization_energy' \
                  ' instead.', DeprecationWarning, stacklevel=2)
    from ase.data.g2 import get_ionization_energy
    return get_ionization_energy(name, vertical)

def get_atomization_energy(name):
    "Deprecated, use ase.data.g2.get_atomization_energy instead."
    import warnings
    warnings.warn('ase.data.molecules.get_atomization_energy is deprecated. '
                  ' Please use from ase.data.g2 import get_atomization_energy' \
                  ' instead.', DeprecationWarning, stacklevel=2)
    from ase.data.g2 import get_atomization_energy
    return get_atomization_energy(name)

def molecule(name, **kwargs):
    "Deprecated."
    import warnings
    warnings.warn('ase.data.molecules.molecule is deprecated. '
                  'Please use:' \
                  ' from ase.structure import molecule' \
                  ' instead.', DeprecationWarning, stacklevel=2)
    from ase.structure import molecule
    return molecule(name, **kwargs)

def latex(name):
    """Convert name to LaTeX"""
    s = '$'
    last = False
    for i in name:
        if i.isalpha():
            if not last:
                s = s + r'\rm{'
                last = True
        elif last:
            s = s + '}'
            last = False
        s = s + i
    if i.isalpha():
        s = s + '}'
    s = s.replace(' ', r'\ ') + '$'
    return s


def rest(name):
    """Convert name to reStructuredText."""
    s = ''
    while name:
        c = name[0]
        if c == '_':
            s += r'\ :sub:`%s`\ ' % name[1]
            name = name[2:]
        elif c == '^':
            s += r'\ :sup:`%s`\ ' % name[1]
            name = name[2:]
        else:
            s += c
            name = name[1:]
    return s
