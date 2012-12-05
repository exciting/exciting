#!/usr/bin/env python

'''
Get the space group for a ListOfAtoms

can be called as a script on a netcdf file

sgroup.py [-f] ncfile

http://cpc.cs.qub.ac.uk/summaries/ADON.html

PROGRAM SUMMARY [Licence| Download | E-mail] adon.tar.gz(69 Kbytes) 
Manuscript Title: Determination of the space group and unit cell for a periodic solid. 
Authors: B.Z. Yanchitsky, A.N. Timoshevskii 
Program title: SGROUP 
Catalogue identifier: ADON 
Journal reference: Comput. Phys. Commun. 139(2001)235 
Programming language: C. 
Computer: Intel/Pentium, Alpha Workstation. 
Operating system: Slackware Linux 4.0, Digitial Unix 4.0D. 
RAM: 1M words 
Word size: 8 
Keywords: Unit cell, Space group, Symmetry operations, Solid state physics, Crystal structure. 
Classification: 7.8. 

'''
import math,os,re,string,tempfile

from Scientific.Geometry import *
from Scientific.IO.FortranFormat import *
from numpy import *

class SGROUP:

    def __init__(self,atoms,outfile=None):
        '''outfile is where the results will be stored if you want
        them. Otherwise they go into a tempfile that is deleted.'''

        id,infile = tempfile.mkstemp()
        
        if outfile is None:
            od,ofile = tempfile.mkstemp()
        else:
            ofile = outfile

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

        f = open(infile,'w')
        f.write('P\n')

        f.write('%1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n' % (a,b,c,
                                                           alpha,beta,gamma))

        f.write('%i\n' % len(atoms))
        
        for i,atom in enumerate(atoms):
            f.write('%1.4f %1.4f %1.4f\n' % tuple(scaledpositions[i]))
            f.write('%s\n\n' % chemicalsymbols[i])
        f.close()

        os.system('sgroup %s %s' % (infile,ofile))
          
        f = open(ofile,'r')
        self.output= f.readlines()
        f.close()

        os.unlink(infile)
        os.close(id)
        
        if outfile is None:
            os.unlink(ofile)
            os.close(od) # you must close the file descriptor or
                         # eventually too many open files will occur
                         # and cause an error when you are processing
                         # many files.

            
    def __str__(self):
        return string.join(self.output)
        
    def get_space_group(self):
        'returns spacegroup number'
        regexp = re.compile('^Number and name of space group:')
        for line in self.output:
            if regexp.search(line):
                line = line[32:]    
                r2 = re.compile('^\d+')
                s = r2.search(line)
                if hasattr(s,'group'):
                    return int(s.group())
                else:
                    return None

    def get_symmetry_operators(self):
        '''
        gets symmetry operators from output

        it looks like this in the output.
        I am not sure what the 4th number is, it is called
        tau in Wien2k. I do not use it or return it here, but
        it is parsed, and could be returned.
        
        Number of symmetry operations: 48
        Operation: 1
        1.0   0.0   0.0  0.000
        0.0   1.0   0.0  0.000
        0.0   0.0   1.0  0.000

        Operation: 2
        -1.0   0.0   0.0  0.000
        0.0  -1.0   0.0  0.000
        0.0   0.0   1.0  0.000

        Operation: 3
        '''
        
        re1 = '^Number of symmetry operations:'
        regexp = re.compile(re1)
        for i,line in enumerate(self.output):
            if regexp.search(line):
                # take integer after the colon
                nsymops = int (string.split(line,':')[-1])
                index = i
                break

        symmetry_operators = []
        taus = []
        for s in range(nsymops):
            temparray=zeros((3,3))
            temptau = [0,0,0]
            if int(string.split(self.output[index+1],':')[-1]) != s+1:
                raise Exception,'this symmetry operator %i does not match index' % s
                               
            x,y,z,tau = [float(var) for var in string.split(self.output[index+2])]
            temparray[0] = [x,y,z]
            temptau[0] = tau
            
            x,y,z,tau = [float(var) for var in string.split(self.output[index+3])]
            temparray[1] = [x,y,z]
            temptau[1] = tau
            
            x,y,z,tau = [float(var) for var in string.split(self.output[index+4])]
            temparray[2] = [x,y,z]
            temptau[2] = tau

            # increase index for next operator
            index += 5

            symmetry_operators.append(temparray)
            taus.append(array(temptau))

        return symmetry_operators,taus
            



if __name__ == '__main__':
    from ase.calculators.jacapo import *
    from optparse import OptionParser

    parser = OptionParser(usage='sgroup.py ncfile',
                      version='0.1')

    parser.add_option('-f',
                      nargs=0,
                      help = 'print full output')

    parser.add_option('-o',
                      nargs=1,
                      help = 'save output in filename')

    options,args = parser.parse_args()

    #print options
    
    for ncfile in args:       

        sg = SGROUP(Jacapo.read_atoms(ncfile),outfile=options.o)

        print sg.get_space_group()
    
        if options.f is not None:
            print sg

