#!/usr/bin/env python
"""
A script which averages a RHO3D.xml or VCL3D.xml file in one direction to make a 1D curve.
User must specify filename and direction on command line.
Depends on lxml
"""

import os
import sys
import numpy as np
import math
import string
import datetime
import time
from lxml import etree as ET

starttime = time.clock() 
print "Starting calculation at",
print time.strftime("%H:%M:%S on %a %d %b %Y")

if len(sys.argv) != 3:
    print "\n** ERROR: Must specify name of file and direction on command line.\n"
    print "** Usage: ", sys.argv[0], " <3D data xml file>  <vacuum direction>"
    sys.exit(0)

if not os.path.isfile(sys.argv[1]):
    print "\n** ERROR: Input file %s was not found." % sys.argv[1]
    sys.exit(0)

# Read information from command line
# First specify location of VCL3D 
VCL3Dfile = sys.argv[1].lstrip()

# Next the direction to make average in 
# input should be x y z, or X Y Z. Default is Z.
allowed = "xyzXYZ"
direction = sys.argv[2].lstrip()
if allowed.find(direction) == -1 or len(direction)!=1 :
    print "** WARNING: The direction was input incorrectly."  
    print "** Setting to z-direction by default."  
if direction.islower():
    direction = direction.upper()
filesuffix = "_%s" % direction

# Open geometry and density class objects
#-----------------------------------------
print "\nReading file: %s" % VCL3Dfile
doc = ET.parse(VCL3Dfile)
gt=[int(n) for n in doc.xpath('//grid/@gridticks')[0].split()]
potl_ = []
row = doc.xpath('//row/row/text()')
for nr in row:
    potl_.append([float(rho) for rho in nr.split()])
potl_ = zip(*[iter(potl_)]*gt[1])

# For VCL3D files we multiply by the volume to get back to eV
Bohr = 0.5291772575
Ha = 27.2113956555
ep = doc.xpath('//axis/@endpointrs')
cell = []
for vec in ep:
    cell.append(vec.split())
cell = np.array(cell).astype(np.float)

# Read in 3D data
#--------------------------------
print "Reading 3D data from file...",
potl_ = np.array(potl_).swapaxes(0,2)

# remove last points from data (to avoid double counting due to PBC)
ngridpts = np.array(potl_.shape)
potl = np.zeros([ngridpts[0]-1,ngridpts[1]-1,ngridpts[2]-1],np.float)
potl = potl_[:ngridpts[0]-1,:ngridpts[1]-1,:ngridpts[2]-1]
del potl_

ngridpts = np.array(potl.shape)
totgridpts = ngridpts.prod()
print "Potential stored on a %dx%dx%d grid" % (ngridpts[0],ngridpts[1],ngridpts[2])
print "Total number of points is %d" % totgridpts

# Find length of lattice vectors
#--------------------------------
latticelength = np.dot(cell, cell.T).diagonal()
latticelength = latticelength**0.5
#print 'latticelength=', latticelength

# End reading file
#------------------------
sys.stdout.flush()
print "done." 

# Perform average
#-----------------
print "Performing average in %s direction" % direction
if direction=="X":
    idir = 0
    a = 1
    b = 2
elif direction=="Y":
    a = 0
    idir = 1
    b = 2
else:
    a = 0
    b = 1
    idir = 2
a = (idir+1)%3
b = (idir+2)%3
# At each point, sum over other two indices
average = np.zeros(ngridpts[idir],np.float)
for ipt in range(ngridpts[idir]):
    if direction=="X":
        average[ipt] = potl[ipt,:,:].sum()
    elif direction=="Y":
        average[ipt] = potl[:,ipt,:].sum()
    else:
        average[ipt] = potl[:,:,ipt].sum()

# Scale by number of grid points in the plane
# The resulting unit will be eV.
average /= ngridpts[a]*ngridpts[b]

# Scale by size of area element in the plane,
# gives unit e/Ang. I.e. integrating the resulting
# RHO file should give the total charge.
# area = np.linalg.det([ (cell[a,a], cell[a,b] ),
#                        (cell[b,a], cell[b,b])])
# dA = area/(ngridpts[a]*ngridpts[b])
# average *= dA

# Print out average
#-------------------
averagefile = VCL3Dfile.split('.')[0] + filesuffix + '.dat'
print "Writing averaged data to file %s..." % averagefile,
sys.stdout.flush()
outputfile = open(averagefile,"w")
outputfile.write("#  Distance (bohr)     Plane averaged potential (eV)\n")

xdiff = latticelength[idir]/float(ngridpts[idir])
for i in range(ngridpts[idir]):
    x = i*xdiff
    outputfile.write("%15.8f %15.8f\n" % (x,average[i]))
outputfile.write("%15.8f %15.8f\n" % (x+xdiff,average[0]))
outputfile.close()
print "done."

endtime = time.clock() 
runtime = endtime-starttime
print "\nEnd of calculation." 
print "Program was running for %.2f seconds." % runtime

sys.exit()
