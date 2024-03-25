#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   sys  import stdin
from   math import sqrt
from   math import factorial
import numpy
import sys
import os

#-------------------------------------------------------------------------------

def fit(order,x,y,npoints,nderiv):
    import math
    b = 99999999
    if (npoints > order): 
       b = numpy.polyfit(x,y,order)[order-nderiv]
       b = math.factorial(nderiv)*b
    return b

def printderiv(etam,bb,nfit,output_file):
    print >> output_file, '%12.8f'%(etam),
    for j in range(0,nfit):
        if (bb[j] != 99999999):
            print >> output_file, '%14.6f'%(bb[j]),
    print >> output_file
    return 

def printnice(etam,bb,nfit,nderiv,orderstep,displpoints):
    invcm2hz = 33.356409
    print
    print "#############################################\n"
    print " Fit data------------------------------------\n"
    print " Maximum value of the displacement  ==>", '%5.3f'%(etam) 
    print " Number of displacement values used ==>", '%5i'%(displpoints),"\n" 
    print " Fit results for the derivative of order", '%4i'%(nderiv),"\n"  
    for j in range(0,nfit):
        order=nderiv+j*orderstep
        if (bb[j] != 99999999):
            print " Polynomial of order", '%2i'%(order), " ==> ", 
            print '%8.2f'%(bb[j]), "[cm-1]"
    print
    for j in range(0,nfit):
        order=nderiv+j*orderstep
        if (bb[j] != 99999999):
            print " Polynomial of order", '%2i'%(order), " ==> ", 
            print '%8.4f'%(bb[j]/invcm2hz), " [THz]"
    print
    return 

def sortdispl(s,e):
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

if (str(os.path.exists('INFO-diamond-phonon'))=='False'): 
    sys.exit("ERROR: file INFO-diamond-phonon not found!\n")
if (str(os.path.exists('force-vs-displacement'))=='False'): 
    sys.exit("ERROR: file force-vs-displacement not found!\n")

input_info = open('INFO-diamond-phonon',"r")
line = input_info.readline()
lattice_parameter = float(input_info.readline().strip().split()[5])
input_info.close()

#-------------------------------------------------------------------------------

maximum_displ = input("\nEnter maximum displacement for the fit >>>> ")
order_of_derivative = input("\nEnter the order of derivative >>>> ")
if (order_of_derivative < 0): 
    sys.exit("ERROR: Order of derivative must be positive!\n")
amass = input("\nEnter atomic mass >>>> ")

#-------------------------------------------------------------------------------

orderstep = 1
factor = -2.
amu = 1.66053886
invcm2hz = 33.356409
bohr_radius = 0.529177
pi = 3.1415926535897931
joule2hartree = 4.3597482
invcm2joule = 5.03411721428
joule2rydberg = joule2hartree/2.

unitconv=invcm2hz/bohr_radius
unitconv=unitconv**2
unitconv=unitconv*factor*joule2hartree/amass/amu/lattice_parameter*10**5

#-------------------------------------------------------------------------------

energy = []
displ = []

#-------------------------------------------------------------------------------

input_energy = open('force-vs-displacement',"r")

while True:
    line = input_energy.readline()
    line = line.strip()
    if len(line) == 0: break
    eta = line.split()[0] 
    ene = line.split()[1] 

    if (abs(float(eta)) <= maximum_displ):  
       energy.append(float(ene))
       displ.append(float(eta))

#-------------------------------------------------------------------------------

displ,energy=sortdispl(displ,energy)

#-------------------------------------------------------------------------------

orderlist=[]

for i in range(0,6):
    dumorder=order_of_derivative+orderstep*i
    orderlist.append(dumorder)   

#-------------------------------------------------------------------------------

output_file = open('check-energy-derivatives',"w")
nmax=len(displ)

while (len(displ) > order_of_derivative and len(displ) > 1): 
    bb = []
    n=len(displ)
    etam=max(displ)
    emin=min(displ)
    etam=max(abs(etam),abs(emin))
    for order in orderlist: 
        cc=fit(order,displ,energy,n,order_of_derivative)
        if (cc != 99999999):
            if (cc >  0): freq=-sqrt(-cc*unitconv)/2./pi
            if (cc <= 0): freq=sqrt( cc*unitconv)/2./pi
            bb.append(freq)
        else:
            bb.append(cc)
    if (n == nmax):
        printnice(etam,bb,6,order_of_derivative,orderstep,nmax)
    printderiv(etam,bb,6,output_file)

    if (abs(displ[0]+etam) < 1.e-7): 
        displ.pop(0)
        energy.pop(0)

    if (abs(displ[len(displ)-1]-etam) < 1.e-7): 
        displ.pop()
        energy.pop()

output_file.close()

output_file = open('order-of-derivative',"w")
print >> output_file, order_of_derivative
output_file.close()
print "#############################################\n"
os.system("PLOT-checkderiv.py")







