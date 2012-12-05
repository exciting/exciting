#!/usr/bin/env python
'''
isotropy 
http://stokes.byu.edu/isolinux.html

http://stokes.byu.edu/iso.tar.gz

You will need to create a directory to unzip the above tarfile in,

cd
mkdir iso
cd iso
wget http://stokes.byu.edu/iso.tar.gz
tar xvzf  iso.tar.gz

#put this in your .cshrc
setenv ISODATA $HOME/iso/
set path=($HOME/iso $path)

'''

import math,os,re,string

from Scientific.Geometry import Vector

class FINDSYM:
    def __init__(self,atoms,outfile=None):
        
        unitcell = atoms.get_cell()
        A = Vector(unitcell[0])
        B = Vector(unitcell[1])
        C = Vector(unitcell[2])

        # lengths of the vectors
        a = A.length()#*angstroms2bohr
        b = B.length()#*angstroms2bohr
        c = C.length()#*angstroms2bohr

        # angles between the vectors
        rad2deg = 360./(2.*math.pi)
        alpha = B.angle(C)*rad2deg
        beta = A.angle(C)*rad2deg
        gamma = A.angle(B)*rad2deg

        scaledpositions = atoms.get_scaled_positions()
        chemicalsymbols = [atom.get_symbol() for atom in atoms]

        input = ''

        input += 'title \n'
        input += '0    tolerance\n'
        input += '2        lattice parameters in lengths and angles\n'
        input += '%1.3f %1.3f %1.3f %1.3f %1.3f %1.3f\n' % (a,b,c,
                                                            alpha,beta,gamma)
        input += '1  3 basis vectors for unit cell\n'

        input += '1.00 0.00 0.00\n'
        input += '0.00 1.00 0.00\n'
        input += '0.00 0.00 1.00\n'

        input += '%i   number of atoms\n' % len(atoms)

        types = ''
        for atom in atoms:
            types += str(atom.get_atomic_number()) + ' '

        input += types + '\n'

        for i,atom in enumerate(atoms):
            input += '%1.3f %1.3f %1.3f\n' % tuple(scaledpositions[i])

        pin,pout = os.popen2('findsym')
        pin.writelines(input)
        pin.close()
        self.output = pout.readlines()
        pout.close()

        if outfile:
            f = open(outfile,'w')
            f.writelines(self.output)
            f.close()

        if os.path.exists('findsym.log'):
            os.remove('findsym.log')

    def __str__(self):
        return string.join(self.output)

    def get_space_group(self):
        regexp = re.compile('^Space Group')

        for line in self.output:
            if regexp.search(line):
                return line


if __name__ == '__main__':
    from ase.calculators.jacapo import *
    from optparse import OptionParser

    parser = OptionParser(usage='findsym.py ncfile',
                      version='0.1')

    parser.add_option('-f',
                      nargs=0,
                      help = 'print full output')

    parser.add_option('-o',
                      nargs=1,
                      help = 'save output in filename')

    options,args = parser.parse_args()
    
    for ncfile in args:       

        sg = FINDSYM(Jacapo.read_atoms(ncfile),outfile=options.o)

        print sg.get_space_group()
    
        if options.f is not None:
            print sg
