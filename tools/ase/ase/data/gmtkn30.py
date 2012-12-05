import os
import pprint
import re
from urllib import urlretrieve
import zipfile
import shutil

import datetime

import numpy as np

from ase.units import Bohr
from ase.atom import Atom
from ase.atoms import Atoms
from ase.data import atomic_numbers, chemical_symbols

# databases from http://toc.uni-muenster.de/GMTKN/GMTKN30/GMTKN30main.html
url_root = 'http://toc.uni-muenster.de/GMTKN/GMTKN30/'
# we may store all downloaded files locally
# (a good idea, but need to ask permission from the authors)
#url_root = './GMTKN30/'
databases = [
    'MB08-165', # 180
    'W4-08', # 111
    'G21IP', # 71
    'G21EA', # 50
    'PA', # 24
    'SIE11', # 29
    'BHPERI', # 61
    'BH76', # 95
    'RSE43', # 88
    'O3ADD6', # 9
    'G2RC', # 47
    'AL2X', # 14
    'NBPRC', # 21
    'ISO34', # 63
    'ISOL22', # 44
    'DC9', # 19
    'DARC', # 22
    'ALK6', # 13
    'BSR36', # 38
    'IDISP', # 13
    'WATER27', # 30
    'S22', # 57
    'ADIM6', # 12
    'RG6', # 11
    'HEAVY28', # 38
    'PCONF', # 11
    'ACONF', # 18
    'SCONF', # 19
    'CYCONF', # 11
    ]

database_files = {}
for db in databases:
    database_files[db] = {
        'structures': 'strucs/' + db + 'structures.zip',
        'ref': db + 'ref.html',
        'module': 'GMTKN30_' + db.replace('-', '_'),
        }
    for xc in ['PBE', 'PBE0', 'SVWN']:
        database_files[db][xc] = 'funcsGMTKN30/' + db + xc + '.html'

def download_file(url, filename, dir='.'):
    # do not mirror subdirectory structure of url
    outfile = os.path.join(dir, os.path.basename(filename))
    urlretrieve(os.path.join(url, filename), outfile)
    return outfile

def read_charge_filter(s):
    try:
        return re.search('\(([-+]\d+)\)', s).group(1)
    except AttributeError:
        return False

def read_charge(filename, dir='.'):
    fh = open(os.path.join(dir, filename), 'rb')
    lines = filter(read_charge_filter, fh.readlines())
    charge = []
    for line in lines:
        sline = line.split()
        charge.append((sline[0],
                       float(re.search('\(([-+]\d+)\)', sline[1]).group(1))))
    fh.close()
    return charge

def read_charges(dirname, dir='.'):
    fullname = os.path.join(dir, dirname)
    for root, dirs, files in os.walk(fullname):
        for file in files:
            if file == 'README': # read charge/number of unpaired electrons file
                return read_charge(file, dir=root)
                break
        else:
            return []

def read_number_of_unpaired_electrons_filter(s):
    try:
        return re.search('\((\d+)\)', s).group(1)
    except AttributeError:
        return False

def read_number_of_unpaired_electrons(filename, dir='.'):
    fh = open(os.path.join(dir, filename), 'rb')
    lines = filter(read_number_of_unpaired_electrons_filter, fh.readlines())
    number_of_unpaired_electrons = []
    for line in lines:
        sline = line.split()
        no_unpaired_electrons = float(re.search('\((\d+)\)', sline[1]).group(1))
        number_of_unpaired_electrons.append((sline[0], no_unpaired_electrons))
    fh.close()
    return number_of_unpaired_electrons

def read_numbers_of_unpaired_electrons(dirname, dir='.'):
    fullname = os.path.join(dir, dirname)
    for root, dirs, files in os.walk(fullname):
        for file in files:
            if file == 'README': # read charge/number of unpaired electrons file
                return read_number_of_unpaired_electrons(file, dir=root)
                break
        else:
            return []

def read_geometry_filter(s):
    return (not s.startswith('$'))

def read_geometry(filename, dir='.'):
    fh = open(os.path.join(dir, filename), 'rb')
    lines = filter(read_geometry_filter, fh.readlines())
    # return geometry in ASE format
    geometry = []
    for line in lines:
        sline = line.split()
        # find chemical symbol (the symbols in the file are lowercase)
        symbol = sline[-1]
        for s in chemical_symbols:
            if symbol == s.lower():
                symbol = s
                break
        geometry.append(Atom(symbol=symbol, position=sline[:-1]))
    fh.close()
    atoms = Atoms(geometry)
    atoms.set_positions(atoms.get_positions()*Bohr) # convert to Angstrom
    return atoms

def read_structures(dirname, dir='.'):
    fullname = os.path.join(dir, dirname)
    geometries = []
    for root, dirs, files in os.walk(fullname):
        for file in files:
            if file != 'README': # skip file
                geometries.append((file, read_geometry(file, dir=root)))
    return geometries

def read_html(filename, dir='.'):
    fh = open(os.path.join(dir, filename), 'rb')
    table = fh.read()
    # extract html table: help from David Landis
    table = table.split('<table')
    table = table[1]
    table = table.split('</table')
    table = table[0]
    # keep field separator tags
    table = table.replace('<tr', ' TTRR <')
    table = table.replace('<td', ' TTDD <')
    # remove the html tags
    #table = re.sub('<[^>]+>', '', table) # wrong
    table = re.sub('<.*?>', '', table)
    # remove end-of-line
    table = re.sub('\n', '', table)
    # split on columns
    table = table.split('TTRR')
    csv = []
    separator = ':' # BHPERI contains chemical names with comas
    ncompounds = 0
    for item in table:
        if item.find('TTDD')!=-1:
            item = item.strip().replace('TTDD', separator)
            # remove the first coma
            item = item[1:]
            litem = []
            for f in item.split(separator):
                fs = f.strip()
                try:
                    v = eval(fs)
                    if fs.isdigit() and str(v) != fs: # e.g. undesirable eval('001') = 1
                        v = fs   
                # string: NameError, .*[+-*], etc: SyntaxError
                except (NameError, SyntaxError):
                    v = fs
                litem.append(v)
            # the number of compounds
            # (exclude reference value and reaction number and divide by 2)
            if ncompounds:
                assert ncompounds == (len(litem)-2)/2, 'Error: number of compounds incorrect for reaction: ' + str(litem[0]) + ' in file: ' + filename
            ncompounds = (len(litem)-2)/2
            # set names of unused compounds to empty string
            for i in range(ncompounds):
                if litem[1+i] == 0: litem[1+i] = ''
            # move the reaction identifier to the end of list
            litem.append(litem.pop(0))
            csv.append(litem)
    fh.close()
    # return the number of compounds per reaction, and the table
    return ncompounds, csv

def table2reference(ncompounds, table):
    # convert from format given by read_html
    reactions = []
    reference = {}
    for r in table:
        reaction_id = r[-1]
        reference[reaction_id] = r[-2]
        stoich = []
        for c in range(ncompounds):
            if r[c] != '': # only defined compounds
                # compound names can have spaces around
                stoich.append((str(r[c]).strip(), r[c+ncompounds]))
        stoich.append(('reaction_id', reaction_id))
        reactions.append(stoich)
    return reference, reactions

def table2results(nsets, table, mode='default'):
    assert mode in ['default', 'D3']
    # convert from format given by read_html
    if mode == 'default':
        index = 0
    else:
        index = nsets
    reference = {}
    for r in table[:-3]: # ignore 3 last rows of statistics
        reaction_id = r[-1]
        if r[index] != '': # only defined compounds
            reference[reaction_id] = r[index]
    return reference

def unzip_file(filename, dir='.'):
    # unzip contents of filename into dir
    fh = open(filename, 'rb')
    z = zipfile.ZipFile(fh)
    if not os.path.isdir(dir):
        os.mkdir(dir)
    for entry in z.namelist():
        # skip spurious zip inside zip files (in HEAVY28)
        if entry.find('.zip') == -1:
            outfile = open(entry, 'wb')
            outfile.write(z.read(entry))
            outfile.close()
    fh.close()

def format_data(database, geometries, no_unpaired_electrons=[], charges=[]):
    "Return data in the custom format.  "
    import numpy as np
    data = {}
    for geometry in geometries:
        system = geometry[0]
        atoms = geometry[1]
        # find the heaviest atom in the system
        heaviest = max([a.number for a in atoms])
        heaviest_index = [a.number for a in atoms].index(heaviest)
        # find number of unpaired electrons
        if system in [s[0] for s in no_unpaired_electrons]:
            magmom = 0
            for s, m in no_unpaired_electrons:
                if system == s:
                    magmom = m
                    break
            magmoms = [0.0 for a in atoms]
            # assume the magnetic moment on the heaviest atom in the system
            # this is incorrect, but is there a better way to set the magnetic moment?
            magmoms[heaviest_index] = float(magmom)
            usemagmoms = np.array(magmoms)
        else:
            usemagmoms = None
        # find charge, put it on the heaviest atom
        if system in [s[0] for s in charges]:
            charge = 0
            for s, c in charges:
                if system == s:
                    charge = c
                    break
            cs = [0.0 for a in atoms]
            cs[heaviest_index] = float(charge)
            usecharges = np.array(cs)
        else:
            usecharges = None
        # populate data
        data[system] = {
            'database': database,
            'name': atoms.get_name(),
            'symbols': ''.join(atoms.get_chemical_symbols()),
            'magmoms': usemagmoms, # None or list
            'charges': usecharges, # None or list
            'positions': atoms.get_positions(),
            }
    return data

def main():
    import os
    if not os.path.isdir('GMTKN30/strucs'):
        os.makedirs('GMTKN30/strucs')
    #for database in ['G2RC', 'WATER27']:
    for database in database_files.keys(): # all databases
        fh = open(database_files[database]['module'].lower() + '.py', 'w')
        fh.write('# Computer generated code! Hands off!\n')
        fh.write('# Generated: ' + str(datetime.date.today()) + '\n')
        fh.write('from numpy import array\n')
        fh.write('data = ')
        data = {} # specification of molecules
        info = {} # reference/calculation info
        # download structures
        file = database_files[database]['structures']
        f = os.path.abspath(download_file(url_root, file, dir='GMTKN30/strucs'))
        fdir = os.path.splitext(os.path.basename(f))[0]
        unzip_file(f, dir=fdir)
        structures = read_structures(fdir)
        no_unpaired_electrons = read_numbers_of_unpaired_electrons(fdir)
        charges = read_charges(fdir)
        # remove temporary directory
        if os.path.isdir(fdir): shutil.rmtree(fdir)
        data = format_data(database, structures, no_unpaired_electrons, charges)
        pprint.pprint(data, stream=fh)
        fh.write('info = ')
        # download reference data
        info = {}
        file = database_files[database]['ref']
        f = download_file(url_root, file, dir='GMTKN30')
        ncompounds, table = read_html(f)
        # transform table into reactions format
        reference, reactions = table2reference(ncompounds, table)
        info['reactions'] = reactions
        info['reaction energy'] = {}
        info['reaction energy']['reference'] = reference
        # download XC results
        for xc in ['PBE', 'PBE0', 'SVWN']:
            file = database_files[database][xc]
            f = download_file(url_root, file, dir='GMTKN30')
            nsets, table = read_html(f)
            # transform table into results format
            reference = table2results(nsets, table)
            info['reaction energy'][xc] = reference
        pprint.pprint(info, stream=fh)
        fh.close()

if __name__ == '__main__':
    main()
