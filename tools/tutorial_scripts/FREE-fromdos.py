#!/usr/bin/python2
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   numpy import *
from   sys   import stdin
import math
import numpy
import sys
import os

#-------------------------------------------------------------------------------
# FUNCTIONS

def f(x,epsmin,epsmax):
    b = 1.
    if (x > epsmin):
        b = 0.
        if (x < epsmax): b = x**2*exp(x)/(exp(x)-1)**2
    return b

#-------------------------------------------------------------------------------

def g(t,freq,epsmin,epsmax,cesp):
    b = epsmax
    if (t > epsmin): b = cesp*freq/t
    return b

#-------------------------------------------------------------------------------

def v(x,epsmin,epsmax):
    b = 1./2.
    if (x > epsmin):
        b = x/2
        if (x < epsmax): b = x/2 + math.log(1.-exp(-x))
    return b

#-------------------------------------------------------------------------------

def leggi(filin,x,y):
    filin.readline()
    while True:
        line = filin.readline().strip() 
        if (len(line) == 0): break
        x.append(float(line.split()[0]))
        y.append(float(line.split()[1]))
    return

#-------------------------------------------------------------------------------

def fvib(t,ome,dos,eunit):
    #    
    hzev = 0.4135717e-14 ; cmev = 0.1239848e-3
    ckev = 8.617e-5      ; cesp = cmev/ckev
    cbec = 1.e6          ; cacm = 1.e-8
    xkb  = 1.3807e-16    ; ab   = 0.529177e0
    hzcm = hzev/cmev     ; hesp = cmev/ckev
    ekjm = 96.4853365    ; hacm = 2.194746e5
    evha = 0.036749309
    ha2ev = 27.211396132
    #
    epsmin = 1.e-15 ; epsmax = 50
    #
    factor = ha2ev
    if ( eunit=='Ha' ): factor = 1.

    heat = 0. 
    fvib = 0. 
    zero = 0. 
    uvib = 0. 
    svib = 0.

    for j in range(len(dos)-1):
        d = ome[j+1]-ome[j]
        x = g(t,ome[j],epsmin,epsmax,hesp)
        fvib = fvib + d*dos[j]*v(x,epsmin,epsmax)*t*ckev*evha
        heat = heat + d*dos[j]*f(x,epsmin,epsmax)
        zero = zero + d*dos[j]*x*t*ckev*evha/2.
        uvib = uvib + d*dos[j]*x/math.tanh(x/2.)*t*ckev*evha/2.
        svib = svib + d*dos[j]*(x/math.tanh(x/2.)/2.-v(x,epsmin,epsmax))*evha*ckev
        tsvib = uvib - fvib
    if (t < 1.e-6): svib = 0. 
    return fvib*factor, uvib*factor, tsvib*factor, svib*factor*1000., heat, zero*factor

#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

if (narg<3): 
    print "\n Incorrect number of arguments. **Usage**:\n\n",
    print " FREE-fromdos.py TMIN TMAX NTSTEPS [EUNIT]\n"
    print " Temperatures should be given in Kelvin\n"
    sys.exit()

tmin = float(sys.argv[1]) ; tmax = float(sys.argv[2]) ; ntpt = int(sys.argv[3])
tstep = (tmax-tmin)/float(ntpt)

eunit ='eV'
if (narg>3): eunit = str(sys.argv[4])
    
#-------------------------------------------------------------------------------

fdos = open("PHDOS.OUT","r") 
ome = [] ; dos = [] ; leggi(fdos,ome,dos)

filename = []
filename.append('F_vib')
filename.append('U_vib')
filename.append('TS_vib')
filename.append('S_vib')
filename.append('C_v')

fout = []
for ifile in filename: fout.append(open(ifile,"w"))

#-------------------------------------------------------------------------------

nstates = 0
for i in range(len(ome)-1):
    delta   = ome[i+1]-ome[i] 
    nstates = nstates+delta*dos[i]
    
hacm = 2.194746e5
for i in range(len(ome)):
    ome[i] = ome[i]*hacm 
    dos[i] = dos[i]/hacm

for i in range(ntpt+1):
    t = tmin+i*tstep
    if (t < 1e-10): t = 1e-10
    elist = fvib(t,ome,dos,eunit)
    for j in range(5):
        print >>fout[j], '%16.8e'%(t), '%16.8e'%(elist[j])
    if ( i==0 ): print "\n Zero-point energy is",'%10.4e'%(elist[5]),"["+eunit+"]\n"

#-------------------------------------------------------------------------------

for i in range(len(filename)): fout[i].close()
 
#-------------------------------------------------------------------------------
    
xlabel = " -lx 'Temperature [K]'"
ylabel = " -ly 'Free energy ["+eunit+"]'"
options = xlabel+ylabel+" -s 0.65 1.3 -mtx 4 -g"
fileslist = "-f F_vib U_vib TS_vib"
os.system("PLOT-files.py "+fileslist+options)
os.system("mv PLOT.png PLOT-free-energies.png")
print " Created PNG output \"PLOT-free-energies.png\""
    
#-------------------------------------------------------------------------------

xlabel = " -lx 'Temperature [K]'"
ylabel = " -ly 'Entropy [m"+eunit+"/K]'"
options = xlabel+ylabel+" -g"
fileslist = "-f S_vib"
os.system("PLOT-files.py "+fileslist+options)
os.system("mv PLOT.png PLOT-entropy.png")
print " Created PNG output \"PLOT-entropy.png\""

#-------------------------------------------------------------------------------

xlabel = " -lx 'Temperature [K]'"
ylabel = " -ly 'Heat capacity [$k_B$]'"
options = xlabel+ylabel+" -g"
fileslist = "-f C_v"
os.system("PLOT-files.py "+fileslist+options)
os.system("mv PLOT.png PLOT-heat-capacity.png")
print " Created PNG output \"PLOT-heat-capacity.png\"\n"

#-------------------------------------------------------------------------------
