from math import sin, cos, pi, sqrt

import numpy as np

from ase.atoms import Atoms, Atom
from ase.units import Bohr, Ry

def read_scf(filename):
    try:
        f = open(filename + '.scf', 'r')
        pip = f.readlines()
        ene = []
        for line in pip:
            if line[0:4] == ':ENE':
                ene.append(float(line[43:59]) * Ry)
        f.close()
        return ene
    except:
        return None

def read_struct(filename, ase = True):
    f = open(filename, 'r')
    pip = f.readlines()
    lattice = pip[1][0:3]
    nat = int(pip[1][27:30])
    cell = np.zeros(6)
    for i in range(6):
        cell[i] = float(pip[3][0 + i * 10:10 + i * 10])
    cell[0:3] = cell[0:3] * Bohr
    if lattice == 'P  ':
        lattice = 'P'
    elif lattice == 'H  ':
        lattice = 'P'
        cell[3:6] = [90.0, 90.0, 120.0]
    elif lattice == 'R  ':
        lattice = 'R'
    elif lattice == 'F  ':
        lattice = 'F'
    elif lattice == 'B  ':
        lattice = 'I'
    elif lattice == 'CXY':
        lattice = 'C'
    elif lattice == 'CXZ':
        lattice = 'B'
    elif lattice == 'CYZ':
        lattice = 'A'
    else:
        print 'TEST needed'
    pos = np.array([])
    atomtype = []
    rmt = []
    neq = np.zeros(nat)
    iline = 4
    indif = 0
    for iat in range(nat):
        indifini = indif
        if len(pos) == 0:
            pos = np.array([[float(pip[iline][12:22]),
                             float(pip[iline][25:35]),
                             float(pip[iline][38:48])]])
        else:
            pos = np.append(pos, np.array([[float(pip[iline][12:22]),
                                            float(pip[iline][25:35]),
                                            float(pip[iline][38:48])]]),
                            axis = 0)
        indif += 1
        iline += 1
        neq[iat] = int(pip[iline][15:17])
        iline += 1
        for ieq in range(1, int(neq[iat])):
            pos = np.append(pos, np.array([[float(pip[iline][12:22]),
                                            float(pip[iline][25:35]),
                                            float(pip[iline][38:48])]]),
                            axis = 0)
            indif += 1
            iline += 1
        for i in range(indif - indifini):
            atomtype.append(pip[iline][0:2].replace(' ', ''))
            rmt.append(float(pip[iline][43:48]))
        iline += 4
    if ase:
        cell2 = coorsys(cell)
        atoms = Atoms(atomtype, pos, pbc = True)
        atoms.set_cell(cell2, scale_atoms = True)
        cell2 = np.dot(c2p(lattice), cell2)
        if lattice == 'R':
            atoms.set_cell(cell2, scale_atoms = True)
        else:
            atoms.set_cell(cell2)
        return atoms
    else:
        return cell, lattice, pos, atomtype, rmt

def write_struct(filename, atoms2 = None, rmt = None, lattice = 'P'):
    atoms=atoms2.copy()
    atoms.set_scaled_positions(atoms.get_scaled_positions())
    f = file(filename, 'w')
    f.write('ASE generated\n')
    nat = len(atoms)
    if rmt == None:
        rmt = [2.0] * nat
    f.write(lattice+'   LATTICE,NONEQUIV.ATOMS:%3i\nMODE OF CALC=RELA\n'%nat)
    cell = atoms.get_cell()
    metT = np.dot(cell, np.transpose(cell))
    cell2 = cellconst(metT)
    cell2[0:3] = cell2[0:3] / Bohr
    f.write(('%10.6f' * 6) % tuple(cell2) + '\n')
    #print atoms.get_positions()[0]
    for ii in range(nat):
        f.write('ATOM %3i: ' % (ii + 1))
        pos = atoms.get_scaled_positions()[ii]
        f.write('X=%10.8f Y=%10.8f Z=%10.8f\n' % tuple(pos))
        f.write('          MULT= 1          ISPLIT= 1\n')
        zz = atoms.get_atomic_numbers()[ii]
        if zz > 71:
            ro = 0.000005 
        elif zz > 36:
            ro = 0.00001
        elif zz > 18:
            ro = 0.00005
        else:
            ro = 0.0001
        f.write('%-10s NPT=%5i  R0=%9.8f RMT=%10.4f   Z:%10.5f\n' %
                (atoms.get_chemical_symbols()[ii], 781, ro, rmt[ii], zz))
        f.write('LOCAL ROT MATRIX:    %9.7f %9.7f %9.7f\n' % (1.0, 0.0, 0.0))
        f.write('                     %9.7f %9.7f %9.7f\n' % (0.0, 1.0, 0.0))
        f.write('                     %9.7f %9.7f %9.7f\n' % (0.0, 0.0, 1.0))
    f.write('   0\n')

def cellconst(metT):
    aa = np.sqrt(metT[0, 0])
    bb = np.sqrt(metT[1, 1])
    cc = np.sqrt(metT[2, 2])
    gamma = np.arccos(metT[0, 1] / (aa * bb)) / np.pi * 180.0
    beta  = np.arccos(metT[0, 2] / (aa * cc)) / np.pi * 180.0
    alpha = np.arccos(metT[1, 2] / (bb * cc)) / np.pi * 180.0
    return np.array([aa, bb, cc, alpha, beta, gamma])

def coorsys(latconst):
    a = latconst[0]
    b = latconst[1]
    c = latconst[2]
    cal = np.cos(latconst[3] * np.pi / 180.0)
    cbe = np.cos(latconst[4] * np.pi / 180.0)
    cga = np.cos(latconst[5] * np.pi / 180.0)
    sal = np.sin(latconst[3] * np.pi / 180.0)
    sbe = np.sin(latconst[4] * np.pi / 180.0)
    sga = np.sin(latconst[5] * np.pi / 180.0)
    return np.array([[a, b * cga, c * cbe],
                     [0, b * sga, c * (cal - cbe * cga) / sga],
                     [0, 0, c * np.sqrt(1 - cal**2 - cbe**2 - cga**2 + 2 * cal * cbe * cga) / sga]]).transpose()

def c2p(lattice):
    # apply as eg. cell2 = np.dot(ct.c2p('F'), cell)
    if lattice == 'P':
        cell = np.eye(3)
    elif lattice == 'F':
        cell = np.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
    elif lattice == 'I':
        cell = np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]])
    elif lattice == 'C':
        cell = np.array([[0.5, 0.5, 0.0], [0.5, -0.5, 0.0], [0.0, 0.0, -1.0]])
    elif lattice == 'R':
        cell = np.array([[2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0], [-1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0], [-1.0 / 3.0, -2.0/3.0, 1.0 / 3.0]])

    else:
        print 'lattice is ' + lattice + '!'
    return cell
