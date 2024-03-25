#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   numpy import *
from   sys   import stdin
from   math  import sqrt
from   math  import factorial
from   math  import exp
from   math  import log
import numpy
import sys
import os
from   scipy import interpolate

#===============================================================================

def leggi(filin,factor):
    f = open(filin,"r") ; x = [] ; y = []
    while True:
        line = f.readline().strip() 
        if (len(line) == 0): break
        x.append(float(line.split()[0])) 
        yy = float(line.split()[1])*factor ; y.append(yy) 
    f.close()
    return x,y

#===============================================================================

def myleggi(filin,factor):
    f = open(filin,"r") ; x = [] ; y = []
    while True:
        line = f.readline().strip() 
        if (len(line) == 0): break
        x.append(float(line.split()[0])) 
        yy = float(line.split()[1])*factor ; y.append(yy) 
    f.close()
    f = open("fvib","w")
    for i in range(len(y)):
        print>>f, x[i],y[i]
    f.close()
    return 
    
#===============================================================================
   
def myfit(i,x,y,order):
    xx = [] ; yy = []
    for j in range(len(x)):
        xx.append(x[j]) ; yy.append(y[j][i])
    pit = numpy.polyfit(xx,yy,order)
    phi = pit[-2] ; chi = pit[-3]*2
    return phi, chi 
    
#===============================================================================
    
def myline(d,n,diag,order,factor):

    x = [] ; y = [] ; t = []  

    for i in range(n):
        ipoint = i-n/2 ; filein = 'f.'+str(ipoint)+'.0'
        if (diag > 0): filein = 'f.'+str(ipoint)+'.'+str(ipoint)
        if (diag < 0): filein = 'f.'+str(ipoint)+'.'+str(-ipoint)
        x.append(ipoint*d) ; t,q = leggi(filein,factor) ; y.append(q)
        
    phi = [] ; chi = []

    for i in range(len(t)):
        phi.append(myfit(i,x,y,order)[0]) ; chi.append(myfit(i,x,y,order)[1])

    return t, phi, chi
    
#===============================================================================

narg  = len(sys.argv)-1

if (narg<3): 
    print "\nIncorrect number of arguments. **Usage**:\n\n",
    print "FREE-fromdos.py TMIN TMAX NTSTEPS\n"
    print "Temperatures should be given in Kelvin\n"
    sys.exit()
    
tmin = float(sys.argv[1]) ; tmax = float(sys.argv[2]) ; ntpt = int(sys.argv[3])
tstep = (tmax-tmin)/float(ntpt) ; temp = []
for i in range(ntpt+1): temp.append(tmin+i*tstep)

#-------------------------------------------------------------------------------
    
print
print '------------------------------------------------------------------------'
print 'Linear thermal-expansion coefficient for cubic systems'
print '------------------------------------------------------------------------'
print
alat = 10.21#input("Enter value for a0 [Bohr] >>>> ")
c11  = 163.1#input("Enter value for c11 [GPa] >>>> ")
c12  = 66.1#input("Enter value for c12 [GPa] >>>> ")
print

#-------------------------------------------------------------------------------
# c11=163.1  ;  c12=66.1  ;  

GPa2au = 3.398827e-5
meV2Ha = 0.036749309e-2

v0  = alat**3/4
c11 = c11*GPa2au
c12 = c12*GPa2au
b0  = (c11+2*c12)/3.

#-------------------------------------------------------------------------------
# READ input files

t = [] ; phi = [] ; chi = []
t = [] ; phd = [] ; chd = []

myleggi("f.0.0",meV2Ha)

t,phi,chi = myline(0.005,5,0,3,meV2Ha)
t,phd,chd = myline(0.005,5,1,3,meV2Ha)

#-------------------------------------------------------------------------------
# WRITE phi

thefiles = open('phis',"w") ; thefile0 = open('phi0',"w")

for i in range(len(t)):
    phi[i] = phi[i]/v0 ; phd[i] = phd[i]/v0
    print >>thefiles, '%16.8e'%(t[i]), '%16.8e'%(phi[i])
    print >>thefile0, '%16.8e'%(t[i]), '%16.8e'%((phi[i]+phd[i]/2)/2)

thefiles.close() ; thefile0.close()

#-------------------------------------------------------------------------------
# WRITE chi

thefiles = open('chis',"w") ; thefile0 = open('chi0',"w")

for i in range(len(t)):
    chi[i] = chi[i]/v0 ; chd[i] = chd[i]/v0
    print >>thefiles, '%16.8e'%(t[i]), '%16.8e'%(chi[i])
    print >>thefile0, '%16.8e'%(t[i]), '%16.8e'%(chd[i])

thefiles.close() ; thefile0.close()

#-------------------------------------------------------------------------------
# WRITE elastic constants

phi1  = [] ; chi11 = [] ; chi12 = []

for i in range(len(t)):
    phi1.append((phi[i]+phd[i]/2)/2)
    chi11.append(chi[i])
    chi12.append(chd[i]/2-chi[i])

fc11 = open('ec-c11',"w") ; fc12 = open('ec-c12',"w") ; fb0 = open('ec-b0',"w")

cc11 = [] ; cc12 = [] ; bulk = []

for i in range(len(t)):
    cc11.append((c11+chi11[i])/GPa2au)
    cc12.append((c12+chi12[i])/GPa2au)
    bulk.append((cc11[i]+2*cc12[i])/3.)
    print >>fc11, '%16.8e'%(t[i]), '%16.8e'%(cc11[i])
    print >>fc12, '%16.8e'%(t[i]), '%16.8e'%(cc12[i])
    print >>fb0,  '%16.8e'%(t[i]), '%16.8e'%(bulk[i])
     
fc11.close() ; fc12.close() ; fb0.close()
   
#-------------------------------------------------------------------------------
# WRITE thermal expansion

thefiles = open('eps-nonlinear',"w") ; thefile0 = open('eps-linear',"w")

epss = []    
eps0 = []

for i in range(len(t)):
    epss.append(-1./3./bulk[i]/GPa2au*phi1[i])  
    eps0.append(-1./3./b0*phi1[i])  
    print >>thefiles, '%16.8e'%(t[i]), '%16.8e'%(epss[i])
    print >>thefile0, '%16.8e'%(t[i]), '%16.8e'%(eps0[i])

thefiles.close() ; thefile0.close()

#-------------------------------------------------------------------------------
# WRITE thermal expansion coefficient

thefiles = open('alp-nonlinear',"w") ; thefile0 = open('alp-linear',"w")

splf = interpolate.splrep(t,epss,s=0)
alps = [] ; alps.append(interpolate.splev(temp,splf,der=1))
splf = interpolate.splrep(t,eps0,s=0)
alp0 = [] ; alp0.append(interpolate.splev(temp,splf,der=1))

for i in range(len(temp)):
    print >>thefiles, '%16.8e'%(temp[i]), '%16.8e'%(alps[0][i])
    print >>thefile0, '%16.8e'%(temp[i]), '%16.8e'%(alp0[0][i])

thefiles.close() ; thefile0.close()

#-------------------------------------------------------------------------------
os.system("$MYSCRIPTS/GNU-two alp-linear alp-nonlinear")
