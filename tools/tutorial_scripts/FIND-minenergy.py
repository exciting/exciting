#!/usr/bin/python
#
################################################################################
#
#_______________________________________________________________________________

from   sys  import stdin
from   math import sqrt
from   math import factorial
import numpy
import sys
import os

#-------------------------------------------------------------------------------

def shell_value(variable,vlist,default):
    v = default
    e = False
    for i in range(len(vlist)):
        if ( vlist[i] == variable ): v = os.environ[variable] ; e = True ; break
    return v, e

#-------------------------------------------------------------------------------

def fitmin(order,x,y,npoints,mxs):
    import math
    b = 99999999
    a = []
    aa= []
    if (npoints > order): 
       aa=numpy.polyfit(x,y,order)
       for i in range(0,order):
           a.append(numpy.polyfit(x,y,order)[i]*(order-i))
       roots=numpy.roots(a)
       for root in roots:
           if (abs(root) < mxs and abs(root.imag) == 0):
               b=root
    return b

def sortstrain(s,e):
    ss=[]
    ee=[]
    ww=[]
    for i in range(len(s)): ww.append(s[i])
    ww.sort()
    for i in range(len(s)):
        ss.append(s[s.index(ww[i])])
        ee.append(e[s.index(ww[i])])
    return ss, ee

#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])
pdefault = shell_value('PDEFAULT',ev_list,"")[1] 

#-------------------------------------------------------------------------------

if (str(os.path.exists('energy-vs-strain'))=='False'): 
    sys.exit("\nERROR: file energy-vs-strain not found!\n")

#-------------------------------------------------------------------------------

maximum_strain = 1
order = 3

if (len(sys.argv) > 1): order = int(sys.argv[1])
if (len(sys.argv) > 2): maximum_strain = int(sys.argv[1])
if (order < 2): sys.exit("\nERROR: polynomial order must be larger than 1!\n")

#-------------------------------------------------------------------------------

energy = []
strain = []

#-------------------------------------------------------------------------------

input_energy = open('energy-vs-strain',"r")

while True:
    line = input_energy.readline()
    line = line.strip()
    if len(line) == 0: break
    eta,ene = line.split() 

    if (abs(float(eta)) <= maximum_strain):  
       energy.append(float(ene))
       strain.append(float(eta))

#-------------------------------------------------------------------------------

strain,energy=sortstrain(strain,energy)

#-------------------------------------------------------------------------------

output_file = open('minimum-energy',"w")

n_ini = len(strain)

while (len(strain) > 1): 
    bb = []
    n=len(strain)
    etam=max(strain)
    emin=min(strain)
    etam=max(abs(etam),abs(emin))
    
    cc=fitmin(order,strain,energy,n,etam)
    if (cc != 99999999):
        print >> output_file, '%12.8f'%(etam), '%12.8f'%(cc.real)

    if (abs(strain[0]+etam) < 1.e-7): 
        strain.pop(0)
        energy.pop(0)

    if (abs(strain[len(strain)-1]-etam) < 1.e-7): 
        strain.pop()
        energy.pop()

output_file.close()

os.system("PLOT-one.py minimum-energy")
#-------------------------------------------------------------------------------







