#!/usr/bin/env python

import os
import sys
import numpy as np
import math
import string

from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import griddata

if len(sys.argv) < 3 or len(sys.argv) > 4 :
    print "\n** ERROR: Must specify name of file and direction on command line.\n"
    print "** Usage: ", sys.argv[0], " <KS band structure> <exciton eigenvector> <keyword>"
    sys.exit(0)
if not os.path.isfile(sys.argv[1]):
    print "\n** ERROR: Input file %s was not found." % sys.argv[1]
    sys.exit(0)

if not os.path.isfile(sys.argv[2]):
    print "\n** ERROR: Input file %s was not found." % sys.argv[2]
    sys.exit(0)
if len(sys.argv) == 4 and str(sys.argv[3]) != 'xas' :
    print "\n** ERROR: Unknown Keyword %s .\n" % sys.argv[3]
    sys.exit(0)
if len(sys.argv) == 4 and str(sys.argv[3]) == 'xas':
	xas=True
else:
	xas=False
#----------------------------------------------
# read Band Structure
#----------------------------------------------
f = open(sys.argv[1],'r')

# dimensions
line = f.readline().split()
nbnd0 = int(line[1])-1 # indexing from 0
nbnd1 = int(line[2])-1 # indexing from 0
nkpt  = int(line[3])

path = np.zeros(nkpt)
vkl  = np.zeros((nkpt,3))

nbnd = nbnd1-nbnd0+1
ene  = np.zeros((nkpt,nbnd))

for ist in xrange(nbnd0,nbnd1+1):
  for ik in xrange(nkpt):
    line = f.readline().split()
    vkl[ik,0]   = float(line[2])
    vkl[ik,1]   = float(line[3])
    vkl[ik,2]   = float(line[4])
    path[ik]    = float(line[5])
    ene[ik,ist-nbnd0] = float(line[6])
  f.readline() # skip line

f.close()

#----------------------------------------------
# read Exciton WF Weights
#----------------------------------------------
f = open(sys.argv[2],'r')

# dimensions
line  = f.readline().split()
npnt_ = int(line[1])
nv    = int(line[2])
v0    = int(line[3])-1 # indexing from 0
nc    = int(line[4])
c0    = int(line[5])-1

if v0 < nbnd0:
  print "ERROR(", sys.argv[0],"): Lowest index mismatch v0 < nbnd0!"
  sys.exit(1)

grid_ = np.zeros((npnt_,3))
data_ = np.zeros((npnt_,nv,nc))

for ip_ in range(npnt_):
  for iv in range(nv):
    for ic in range(nc):
      line            = f.readline().split()
      grid_[ip_,0]     = float(line[1])
      grid_[ip_,1]     = float(line[2])
      grid_[ip_,2]     = float(line[3])
      data_[ip_,iv,ic] = float(line[6])
f.close()

# single band contributions
if xas:
	tmp = np.zeros((npnt_,nc))
	for ip_ in range(npnt_):
  		for ic in xrange(nc):
			tmp[ip_,ic] = np.sum(data_[ip_,:,ic])

else: 
	tmp = np.zeros((npnt_,nv+nc))
	for ip_ in range(npnt_):
  		for iv in xrange(nv):
			tmp[ip_,iv] = np.sum(data_[ip_,iv,:])
  		for ic in xrange(nc):
			tmp[ip_,nv+ic] = np.sum(data_[ip_,:,ic])

del data_
data_ = tmp
del tmp

#--------------------------------------------------------
# create a supercell (to cover the interpolation volume)
#--------------------------------------------------------
npnt = 3*3*3*npnt_
grid = np.zeros((npnt,3))
if xas:
	data = np.zeros((npnt,nc))
else:
	data = np.zeros((npnt,nv+nc))
ip = 0
for i in xrange(-1,2):
  for j in xrange(-1,2):
    for k in xrange(-1,2):
      for ip_ in xrange(npnt_):
        grid[ip,0] = grid_[ip_,0]+float(i)
        grid[ip,1] = grid_[ip_,1]+float(j)
        grid[ip,2] = grid_[ip_,2]+float(k)
        data[ip,:] = data_[ip_,:]
        ip += 1

#----------------------------------------------
# Output
#----------------------------------------------
ha2ev = 27.2114
if (xas):
	for ist in xrange(nbnd0,nbnd1+1):
		for ik in xrange(nkpt):
			if (ist>=c0) and (ist<=c0+nc-1):
				xyz = tuple(vkl[ik,:])
				val = griddata(grid, data[:,ist-c0], xyz, method='linear')
			else:
				val = 0.0
			print path[ik], ene[ik,ist-nbnd0]*ha2ev, val
		print "\n"
else:	
	for ist in xrange(nbnd0,nbnd1+1):
		for ik in xrange(nkpt):
			if (ist>=v0) and (ist<=c0+nc-1):
				xyz = tuple(vkl[ik,:])
				val = griddata(grid, data[:,ist-v0], xyz, method='linear')
			else:
				val = 0.0
			print path[ik], ene[ik,ist-nbnd0]*ha2ev, val
		print "\n"

#
# Gnuplot command
#
# p "<output>" u 1:2 lt 3 w l, \
#   "<output>" u 1:2:($3*8) with points lt 1 pt 6 ps variable lc rgb "red"
#
