#!/usr/bin/python
#
################################################################################
#
# basic script to visualize phonon modes
# produces an XYZ and an AXSF file containing a number of
# of configurations to be specified via the
# input of nsteps
# it reads the files PHONON.OUT and input.xml, and the
# species file specified therein
#
# a supercell of n1*n2*n3 is plotted
# change the default of 6*6*6 by specifying three integer values on the command line
#
# the files can be visualized by e.g. XCrysDen or VMD
#
# Revision history:
# Dec 2012, script created (STK)
# Oct 2012, bugfix for crystal/@scale not specified, and changes for molecules (STK)
# Oct 2013, changed output format, species file read in $EXCITINGROOT/species (PP)
# Jul 2014, attribute stretch is now considered, input of displacement scaling from command line (STK)
#
#-------------------------------------------------------------------------------

import sys
import os
from lxml import etree
import numpy
from pylab import *
from scipy import interpolate

#-------------------------------------------------------------------------------
# read environment variable

def shell_value(variable,vlist,default):
    v = default
    e = False
    for i in range(len(vlist)):
        if ( vlist[i] == variable ): v = os.environ[variable] ; e = True ; break
    return v, e

#-------------------------------------------------------------------------------
# parameters

bohr2ang = 0.52917721092
ha2rcm = 2.194746313708e05
twopi = 2.0 * 3.14159265359

#-------------------------------------------------------------------------------
# supercell dimensions and scaling, defaults

n1 = 6  ;  n2 = 6  ;  n3 = 6
scaling = 10.

# optional input for non-default parameters, from command line

narg = len(sys.argv)-1

if (narg == 1):
    scaling = float(sys.argv[1])

if (narg == 3):
    n1 = int(sys.argv[1])  ;  n2 = int(sys.argv[2])  ;  n3 = int(sys.argv[3])

if (narg == 4):
    n1 = int(sys.argv[1])  ;  n2 = int(sys.argv[2])  ;  n3 = int(sys.argv[3])  ;  scaling = float(sys.argv[4])

imt = '%5i'    
print "\n############################################"

#-------------------------------------------------------------------------------
# reading number of steps

nsteps = 1
print "\n--------------------------------------------"

#-------------------------------------------------------------------------------
# reading input.xml

bvec = [[0.0 for x in xrange(3)] for x in xrange(3)]
if (str(os.path.exists('input.xml'))=='False'): 
    sys.exit("\nERROR: file input.xml not found!\n")
tree = etree.parse('input.xml')

if (tree.xpath('/input/structure/crystal/@scale') == []):
    alat = 1.0
else:
    alat = float(tree.xpath('/input/structure/crystal/@scale')[0])

str_stretch = tree.xpath('/input/structure/crystal/@stretch')
if (str_stretch ==[]):
    stretch = [1.,1.,1.]
else: stretch=numpy.array(map(float,str_stretch[0].split()))

blabel = "/input/structure/crystal/basevect"
bvec[0][0],bvec[1][0],bvec[2][0] = [float(x) for x in tree.xpath(blabel+'[1]')[0].text.split()]
bvec[0][1],bvec[1][1],bvec[2][1] = [float(x) for x in tree.xpath(blabel+'[2]')[0].text.split()]
bvec[0][2],bvec[1][2],bvec[2][2] = [float(x) for x in tree.xpath(blabel+'[3]')[0].text.split()]

if (tree.xpath('/input/structure/@cartesian') == []):
    cartesian = False
else:
    if (tree.xpath('/input/structure/@cartesian')[0] == 'true'):
        cartesian = True
    else:
        cartesian = False
if (cartesian):
    print 'Atomic positions are in cartesian coordinates\n'
   
#-------------------------------------------------------------------------------

for i in range(3):
    for j in range(3):
        bvec[i][j] = bvec[i][j]*alat*bohr2ang*stretch[j]

nspec = int(tree.xpath('count(/input/structure/species)'))
natom = [float(x) for x in xrange(nspec)]
natomax = 0
natmtot = 0
for isp in range(nspec):
    s = 'count(/input/structure/species[' + str(isp+1) + ']/atom)'
    natom[isp] = int(tree.xpath(s))
    if natom[isp] > natomax: natomax = natom[isp]
    natmtot = natmtot + natom[isp]

#-------------------------------------------------------------------------------

atcoord = [[[0.0 for x in xrange(3)] for x in xrange(natomax)] for x in xrange(nspec)]
on = [0 for x in xrange(nspec)]
mass = [0.0 for x in xrange(nspec)]
symb = ['' for x in xrange(nspec)]

#-------------------------------------------------------------------------------

ev_list = os.environ.keys()
excitingroot = shell_value('EXCITINGROOT',ev_list,"")[0]
speciespath = excitingroot+"/species"
if (excitingroot==''): 
    sys.exit("\nERROR: environment variable EXCITINGROOT not defined!\n")

for isp in range(nspec):
    for iat in range(natom[isp]):
        s = '/input/structure/species[' + str(isp+1) + ']/atom[' + str(iat+1) + ']/@coord'
        atcoord[isp][iat][0],atcoord[isp][iat][1],atcoord[isp][iat][2] = \
                          [float(x) for x in tree.xpath(s)[0].split()]
    s = '/input/structure/species[' + str(isp+1) + ']/@speciesfile'
    speciesfile = tree.xpath(s)[0]
    s = speciespath + '/' + speciesfile
    stree = etree.parse(s)
    on[isp] = -float(stree.xpath('/spdb/sp/@z')[0])
    mass[isp] = float(stree.xpath('/spdb/sp/@mass')[0])
    symb[isp] = stree.xpath('/spdb/sp/@chemicalSymbol')[0]
    #print 'species ',isp+1,' contains ',natom[isp],' ',symb[isp],' atoms'

print 'Total number of atoms    :', (imt%natmtot)

#-------------------------------------------------------------------------------

qvec = [[0.0 for x in xrange(3)] for x in xrange(100)]
evr = [[[[[0.0 for x in xrange(3)] for x in xrange(natomax)] for x in xrange(nspec)] for x in xrange(3*natmtot)] for x in xrange(100)]
evi = [[[[[0.0 for x in xrange(3)] for x in xrange(natomax)] for x in xrange(nspec)] for x in xrange(3*natmtot)] for x in xrange(100)]
freq = [[0.0 for x in xrange(3*natmtot)] for x in xrange(100)]
dispr = [[[0.0 for x in xrange(3)] for x in xrange(natomax)] for x in xrange(nspec)]
dispi = [[[0.0 for x in xrange(3)] for x in xrange(natomax)] for x in xrange(nspec)]
datcoord = [[[0.0 for x in xrange(3)] for x in xrange(natomax)] for x in xrange(nspec)]
n = [0 for x in xrange(3)]

#-------------------------------------------------------------------------------

if (str(os.path.exists('PHONON.OUT'))=='False'): 
    sys.exit("\nERROR: file PHONON.OUT not found!\n")

f = open('PHONON.OUT','r')
linelist = f.readlines()
filelen = len(linelist)
iline = 0
while iline < filelen:
# first line is empty
    iline += 1
# read number of q-point and q vector
    iqpt = int(linelist[iline].split()[0])
    qvec[iqpt-1][0] = float(linelist[iline].split()[1])
    qvec[iqpt-1][1] = float(linelist[iline].split()[2])
    qvec[iqpt-1][2] = float(linelist[iline].split()[3])
# next line is empty
    iline += 2
    for i in range(3*natmtot):
# read mode number and frequency
        imode = int(linelist[iline].split()[0])
        freq[iqpt-1][imode-1] = float(linelist[iline].split()[1])
        iline += 1
# read eigen vectors
        for isp in range(nspec):
            for iat in range(natom[isp]):
                for ipol in range(3):
                    evr[iqpt-1][imode-1][isp][iat][ipol] = float(linelist[iline].split()[3])
                    evi[iqpt-1][imode-1][isp][iat][ipol] = float(linelist[iline].split()[4])
                    iline += 1
# two empty lines follow
        iline +=1
nqpt = iqpt

#-------------------------------------------------------------------------------

print "Total number of q-points :", (imt%nqpt)
print "--------------------------------------------\n"
fmt = '%10.4f'
imt = '%3i'
for iqpt in range(nqpt):
    print 'q-point',(imt%(iqpt+1)),' is ',qvec[iqpt][0],qvec[iqpt][1],qvec[iqpt][2]
    for imode in range(3*natmtot):
        print '   mode',(imt%(imode+1)),'     frequency:',
        print (fmt%(freq[iqpt][imode]*ha2rcm)),' cm-1'
    print
      
print "############################################\n"

#-------------------------------------------------------------------------------
