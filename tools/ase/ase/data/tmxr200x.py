import os
import pprint
import re
from urllib import urlretrieve
import zipfile

import datetime

import numpy as np

from ase.utils import popen3

import ase.io
from ase.atom import Atom
from ase.atoms import Atoms
from ase.data import atomic_numbers, chemical_symbols
from ase.data import ground_state_magnetic_moments

# Transition Metals First-row (TM1R): 10.1021/ct6001187 # 32 compounds
# Transition Metals Second-row (TM2R): 10.1021/ct700178y # 19 compounds
# Transition Metals Third-row (TM3R): 10.1021/ct800172j # 25 compounds

#http://pubs.acs.org/doi/suppl/10.1021/ct6001187/suppl_file/ct6001187-file002.pdf
#http://pubs.acs.org/doi/suppl/10.1021/ct700178y/suppl_file/ct700178y-file002.pdf
#http://pubs.acs.org/doi/suppl/10.1021/ct800172j/suppl_file/ct800172j_si_001.pdf

url_root = 'http://pubs.acs.org/doi/suppl/'
journal = '10.1021'
database_files = {
    'TM1R2006': {'doi': journal + '/ct6001187', 'module': 'TMXR200X_TM1R2006'},
    'TM2R2007': {'doi': journal + '/ct700178y', 'module': 'TMXR200X_TM2R2007'},
    'TM3R2008': {'doi': journal + '/ct800172j', 'module': 'TMXR200X_TM3R2008'},
    }

database_files['TM1R2006']['pdf'] = database_files['TM1R2006']['doi'] + '/suppl_file/ct6001187-file002.pdf'

database_files['TM2R2007']['pdf'] = database_files['TM2R2007']['doi'] + '/suppl_file/ct700178y-file002.pdf'

database_files['TM3R2008']['pdf'] = database_files['TM3R2008']['doi'] + '/suppl_file/ct800172j_si_001.pdf'

def download_file(url, filename, dir='.'):
    # do not mirror subdirectory structure of url
    outfile = os.path.join(dir, os.path.basename(filename))
    if 0: # fails, use files from disk
        urlretrieve(os.path.join(url, filename), outfile)
    return outfile

def read_geometries(filename, dir='.'):
    txt = os.path.join(dir, filename)
    fh = open(txt, 'rb')
    table = fh.read()
    firstsplit = '(in xyz format):' # TM1R2006 and TM2R2007
    dataformat = 'xyz'
    if table.find('(Gaussian archive entries):') != -1:
        firstsplit = '(Gaussian archive entries):' # TM3R2008
        dataformat = 'gaussian'
    table = table.split(firstsplit)
    table = table[1]
    # remove one or two digit numbers (page numbers/numbers of atoms in xyz format)
    table = re.sub('\n\d\d\n', '\n', table)
    table = re.sub('\n\d\n', '\n', table)
    # remove S + two digit numbers (page numbers)
    table = re.sub('\nS\d\d\n', '\n', table)
    # remove S + one digit (page numbers)
    table = re.sub('\nS\d\n', '\n', table)
    # remove empty lines
    # http://stackoverflow.com/questions/1140958/whats-a-quick-one-liner-to-remove-empty-lines-from-a-python-string
    table = os.linesep.join([s for s in table.splitlines() if s])
    geometries = []
    if dataformat == 'xyz':
        # split on new lines
        table = table.split('\n')
        # mark compound names with ':' tags
        for n, line in enumerate(table):
            if not (line.find('.') != -1):
                # remove method/basis set information
                table[n] = table[n].replace(' BP86/qzvp', '')
                table[n] = ':' + table[n] + ':'
        table = '\n'.join([s for s in table])
        # split into compounds
        # http://simonwillison.net/2003/Oct/26/reSplit/
        # http://stackoverflow.com/questions/647655/python-regex-split-and-special-character
        table = re.compile('(:.*:)').split(table)
        # remove empty elements
        table = [l.strip() for l in table]
        table = [l for l in table if len(l) > 1]
        # extract compounds
        for n in range(0, len(table), 2):
            compound = table[n].replace(':', '').replace(' ', '_')
            geometry = []
            for atom in table[n+1].split('\n'):
                geometry.append(Atom(symbol=atom.split()[0], position=atom.split()[1:]))
            atoms = Atoms(geometry)
            # set the charge and magnetic moment on the heaviest atom (better ideas?)
            heaviest = max([a.get_atomic_number() for a in atoms])
            heaviest_index = [a.get_atomic_number() for a in atoms].index(heaviest)
            charge = 0.0
            if abs(charge) > 0.0:
                charges = [0.0 for a in atoms]
                charges[heaviest_index] = charge
                atoms.set_charges(charges)
            if compound in [ # see corresponding articles
                'Ti(BH4)3',  # TM1R2006
                'V(NMe2)4',  # TM1R2006
                'Cu(acac)2',  # TM1R2006
                'Nb(Cp)(C7H7)_Cs', # TM2R2007
                'CdMe_C3v', # TM2R2007
                ]:
                multiplicity = 2.0
            else:
                multiplicity = 1.0
            if multiplicity > 1.0:
                magmoms = [0.0 for a in atoms]
                magmoms[heaviest_index] = multiplicity - 1
                atoms.set_initial_magnetic_moments(magmoms)
            geometries.append((compound, atoms))
    elif dataformat == 'gaussian':
        # remove new lines
        table = table.replace('\n', '')
        # fix: MeHg(Cl) written as MeHg(CN)
        table = table.replace(
            'MeHg(CN), qzvp (SDD/def-qzvp for metal)\\\\0,1\\Hg,0.,0.,0.1975732257',
            'MeHg(Cl), qzvp (SDD/def-qzvp for metal)\\\\0,1\\Hg,0.,0.,0.1975732257')
        # split on compound end marks
        table = table.split('\\\@')
        # remove empty elements
        table = [l.strip() for l in table]
        table = [l for l in table if len(l) > 1]
        # extract compounds
        for n, line in enumerate(table):
            # split on gaussian separator '\\'
            entries = line.split('\\\\')
            compound = entries[2].split(',')[0].split(' ')[0]
            # charge and multiplicity from gaussian archive
            charge, multiplicity = entries[3].split('\\')[0].split(',')
            charge = float(charge)
            multiplicity = float(multiplicity)
            if compound in ['Au(Me)PMe3']: # in gzmat format!
                # check openbabel version (babel >= 2.2 needed)
                cmd = popen3('babel -V')[1]
                output = cmd.read().strip()
                cmd.close()
                v1, v2, v3 = output.split()[2].split('.')
                v1, v2, v3 = int(v1), int(v2), int(v3)
                if not (v1 > 2 or ((v1 == 2) and (v2 >= 2))):
                    print compound + ': skipped - version of babel does not support gzmat format'
                    continue # this one is given in z-matrix format
                finame = compound.replace('(', '').replace(')', '') + '.orig'
                foname = finame.split('.')[0] + '.xyz'
                fi = open(finame, 'w')
                fo = open(foname, 'w')
                if 1: # how to extract zmat by hand
                    zmat = ['#'] # must start with gaussian input start
                    zmat.extend('@') # separated by newline
                    zmat.extend([compound])
                    zmat.extend('@') # separated by newline
                    zmat.extend([str(int(charge)) + ' ' + str(int(multiplicity))])
                    zmat.extend(entries[3].replace(',', ' ').split('\\')[1:])
                    zmat.extend('@') # atom and variable definitions separated by newline
                    zmat.extend(entries[4].split('\\'))
                    zmat.extend('@') # end with newline
                    for l in zmat:
                        fi.write(l.replace('@', '').replace('=', ' ') + '\n')
                    fi.close()
                if 0:
                    # or use the whole gausian archive entry
                    entries = ''.join(entries)
                    fi.write(entries)
                # convert gzmat into xyz using openbabel (babel >= 2.2 needed)
                cmd = popen3('babel -i gzmat ' + finame + ' -o xyz ' + foname)[2]
                error = cmd.read().strip()
                cmd.close()
                fo.close()
                if not (error.find('0 molecules') != -1):
                    atoms = ase.io.read(foname)
                else:
                    print compound + ': babel conversion failed'
                    continue # conversion failed
            else:
                positions = entries[3].replace(',', ' ').split('\\')[1:]
                geometry = []
                for k, atom in enumerate(positions):
                    geometry.append(Atom(symbol=atom.split()[0],
                                         position=[float(p) for p in atom.split()[1:]]))
                atoms = Atoms(geometry)
            #
            # set the charge and magnetic moment on the heaviest atom (better ideas?)
            heaviest = max([a.get_atomic_number() for a in atoms])
            heaviest_index = [a.get_atomic_number() for a in atoms].index(heaviest)
            if abs(charge) > 0.0:
                charges = [0.0 for a in atoms]
                charges[heaviest_index] = charge
                atoms.set_charges(charges)
            if multiplicity > 1.0:
                magmoms = [0.0 for a in atoms]
                magmoms[heaviest_index] = multiplicity - 1
                atoms.set_initial_magnetic_moments(magmoms)
            geometries.append((compound, atoms))
    return geometries

def pdftotext(filename):
    os.system('pdftotext -raw -nopgbrk '+ filename)
    return os.path.splitext(filename)[0] + '.txt'

from ase.data.gmtkn30 import format_data

def main():
    if not os.path.isdir('TMXR200X'):
        os.makedirs('TMXR200X')
    #for database in ['TM1R2006']:
    for database in database_files.keys():
        fh = open(database_files[database]['module'].lower() + '.py', 'w')
        fh.write('# Computer generated code! Hands off!\n')
        fh.write('# Generated: ' + str(datetime.date.today()) + '\n')
        fh.write('from numpy import array\n')
        fh.write('data = ')
        data = {} # specification of molecules
        info = {} # reference/calculation info
        # download structures
        file = database_files[database]['pdf']
        f = os.path.abspath(download_file(url_root, file, dir='TMXR200X'))
        f = pdftotext(f)
        geometries = read_geometries(f)
        # set number of unpaired electrons and charges
        no_unpaired_electrons = []
        charges = []
        for a in geometries:
            magmom = sum(a[1].get_initial_magnetic_moments())
            if magmom > 0.0:
                no_unpaired_electrons.append((a[0], magmom))
            charge = sum(a[1].get_charges())
            if abs(charge) > 0.0:
                charges.append((a[0], charge))
        data = format_data(database, geometries, no_unpaired_electrons, charges)
        # all constituent atoms
        atoms = []
        for formula, geometry in geometries:
            atoms.extend(list(set(geometry.get_chemical_symbols())))
        atoms=list(set(atoms))
        atoms.sort()
        for atom in atoms:
            magmom=ground_state_magnetic_moments[atomic_numbers[atom]]
            data[atom] = {
                'database': database,
                'name': atom,
                'symbols': atom,
                'magmoms': [magmom], # None or list
                'charges': None, # None or list
                'positions': np.array([[0.0]*3]),
                }
            Atom(atom, magmom=magmom)
        pprint.pprint(data, stream=fh)
        fh.close()

if __name__ == '__main__':
    main()
