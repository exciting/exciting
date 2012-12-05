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

#-------------------------------------------------------------------------------
# FUNCTIONS

def f(x,epsmin,epsmax):
    import math
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
    import math
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

def fvib(t,ome,dos,hesp,ckev,evha):
    import math
    epsmin = 1.e-15 ; epsmax = 50
    c = 0. ; e = 0. ; z = 0. ; i = 0. ; s = 0.
    for j in range(len(dos)-1):
        d = ome[j+1]-ome[j]
        x = g(t,ome[j],epsmin,epsmax,hesp)
        e = e+d*dos[j]*v(x,epsmin,epsmax)*t*ckev*evha
        c = c+d*dos[j]*f(x,epsmin,epsmax)
        z = z+d*dos[j]*x*t*ckev*evha/2.
        i = i+d*dos[j]*x/math.tanh(x/2.)*t*ckev*evha/2.
        s = s+d*dos[j]*(x/math.tanh(x/2.)/2.-v(x,epsmin,epsmax))*evha*ckev
    if (t < 1.e-6): s = 0. 
    return e,c,z,s,i

#-------------------------------------------------------------------------------
# GENERAL DATA

hzev = 0.4135717e-14 ; cmev = 0.1239848e-3
ckev = 8.617e-5      ; cesp = cmev/ckev
cbec = 1.e6          ; cacm = 1.e-8
xkb  = 1.3807e-16    ; ab   = 0.529177e0
hzcm = hzev/cmev     ; hesp = cmev/ckev
ekjm = 96.4853365    ; hacm = 2.194746e5
evha = 0.036749309

#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

if (narg<3): 
    print "\nIncorrect number of arguments. **Usage**:\n\n",
    print "FREE-fromdos.py TMIN TMAX NTSTEPS\n"
    print "Temperatures should be given in Kelvin\n"
    sys.exit()

tmin = float(sys.argv[1]) ; tmax = float(sys.argv[2]) ; ntpt = int(sys.argv[3])
tstep = (tmax-tmin)/float(ntpt)

#-------------------------------------------------------------------------------

fdos = open("PHDOS.OUT","r") ; ome = [] ; dos = [] ; leggi(fdos,ome,dos)

print
print "-----------------------------------"
print " What do you want to calculate?"
print "-----------------------------------"
print "  0 =>  Zero-point energy"
print "-----------------------------------"
print "  1 =>  Vibrational internal energy"
print "  2 =>  Vibrational free energy"
print "  3 =>  Vibrational entropy"
print "  4 =>  Heat capacity"
print "-----------------------------------"

task = raw_input("\nEnter task code >>>> ")

print
if (task != "1" and task != "2" and\
    task != "3" and task != "4" and task != "0"): 
    sys.exit("ERROR: Task code is out of range [0-4]!\n")
    
task=int(task)

if (task == 1): filename = 'vibrational-internal-energy'
if (task == 2): filename = 'vibrational-free-energy'
if (task == 3): filename = 'vibrational-entropy'
if (task == 4): filename = 'heat-capacity'

if (str(task) != "0"): fout = open(filename,"w")

#-------------------------------------------------------------------------------

nstates = 0
for i in range(len(ome)-1):
    delta   = ome[i+1]-ome[i] 
    nstates = nstates+delta*dos[i]

for i in range(len(ome)):
    ome[i] = ome[i]*hacm 
    dos[i] = dos[i]/hacm

for i in range(ntpt+1):
    t = tmin+i*tstep
    if (t < 1e-10): t = 1e-10
    ven,hca,zpe,ent,ien = fvib(t,ome,dos,hesp,ckev,evha)
    if (task == 1): print >>fout, '%16.8e'%(t), '%16.8e'%(ien)
    if (task == 2): print >>fout, '%16.8e'%(t), '%16.8e'%(ven)
    if (task == 3): print >>fout, '%16.8e'%(t), '%16.8e'%(ent)
    if (task == 4): print >>fout, '%16.8e'%(t), '%16.8e'%(hca)
    
if (task == 0): print "Zero-point energy is ",'%14.8e'%(zpe)," [Ha]\n"

#-------------------------------------------------------------------------------

if (str(task) != "0"):
    fout.close()
    style    = " o-"
    if (ntpt > 40): style = "-"
    xlabel   = " \"Temperature [K]\""
    if (task == 1): ylabel = " \"Energy [Ha]\"" 
    if (task == 2): ylabel = " \"Free energy [Ha]\"" 
    if (task == 3): ylabel = " \"Entropy [Ha/K]\"" 
    if (task == 4): ylabel = " \"Heat capacity [kB]\"" 
    options = style+xlabel+ylabel
    os.system("PLOT-plot.py "+filename+" "+options)
    print "Created PostScript output \"PLOT.ps\" from file \""+filename+"\"\n"

#-------------------------------------------------------------------------------

