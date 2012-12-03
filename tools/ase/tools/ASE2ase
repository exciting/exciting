#!/usr/bin/env python

from optparse import OptionParser

description = """Convert ASE2 script FILEs to ase3.  FILEs will be
modified in-place to be compatible with ase3.  Original files are
backed up."""

p = OptionParser(usage='%prog FILE...', description=description)

opts, args = p.parse_args()

def convert(filename):
    lines = open(filename).readlines()
    t1 = ''.join(lines)

    first = True
    for i in range(len(lines)):
        line = lines[i]
        if line.startswith('from ASE'):
            if first:
                lines[i] = 'from ase.all import *\n'
                first = False
            else:
                lines[i] = ''

    t = ''.join(lines)

    for old, new in [('GetCartesianPositions', 'get_positions'),
                     ('SetCartesianPositions', 'set_positions'),
                     ('GetPotentialEnergy', 'get_potential_energy'),
                     ('SetCalculator', 'set_calculator'),
                     ('GetScaledPositions', 'get_scaled_positions'),
                     ('SetScaledPositions', 'set_scaled_positions'),
                     ('SetUnitCell', 'set_cell'),
                     ('GetUnitCell', 'get_cell'),
                     ('GetBoundaryConditions', 'get_pbc'),
                     ('GetCartesianForces', 'get_forces'),
                     ('GetCartesianVelocities', 'get_velocities'),
                     ('SetCartesianVelocities', 'set_velocities'),
                     ('GetCartesianMomenta', 'get_momenta'),
                     ('SetCartesianMomenta', 'set_momenta'),
                     ('ListOfAtoms', 'Atoms'),
                     ('periodic', 'pbc'),
                     ('pbcity', 'periodicity'),
                     ('.Converge(', '.run('),
                     ('Repeat', 'repeat'),
                     ('Numeric', 'numpy'),
                     ('numpyal', 'Numerical'),
                     ('GetAtomicNumber()', 'number'),
                     ('GetChemicalSymbol()', 'symbol'),
                     ('GetCartesianPosition()', 'position'),
                     ('GetTag()', 'tag'),
                     ('GetCharge()', 'charge'),
                     ('GetMass()', 'mass'),
                     ('GetCartesianMomentum()', 'momentum'),
                     ('GetMagneticMoment()', 'magmom'),
                     ]:
        t = t.replace(old, new)

    t2 = ''
    while 1:
        i = t.find('.')
        i2 = t.find('def ')
        if 0 <= i < i2:
            n = 1
        elif i2 != -1:
            n = 4
            i = i2
        else:
            break
        t2 += t[:i + n]
        t = t[i + n:]
        if t[0].isupper() and t[1].islower():
            j = t.find('(')
            if j != -1 and t[2: j].isalpha():
                for k in range(j):
                    if t[k].isupper() and k > 0:
                        t2 += '_'
                    t2 += t[k].lower()
                t = t[j:]

    t2 += t

    if t2 != t1:
        print filename, len(t1) - len(t2)
        open(filename + '.bak', 'w').write(t1)
        open(filename, 'w').write(t2)

for filename in args:
    convert(filename)
