#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

import os
import sys  
import math 
import numpy
from   scipy import *
from   scipy.optimize import leastsq

#-------------------------------------------------------------------------------

def inifit(x,y):
    p=[] ; a,b,c = polyfit(x,y,2)
    p.append(-b/(2*a))               # v0
    p.append(a*p[0]**2 + b*p[0] + c) # e0
    p.append((2*a*p[0]))             # b0
    p.append(2.)                     # bp
    return p

def ameos(v,p):
    v0 = p[0] ; e0 = p[1] ; b0 = p[2] ; bp = p[3]  
    vv = (v0/v)**(2./3)
    ff = vv - 1.
    ee = e0 + 9.*b0*v0/16. * ( ff**3*bp + ( 6.-4.*vv )*ff**2 )
    return ee   

def bmeos(v,p):
    v0 = p[0] ; e0 = p[1] ; b0 = p[2] ; bp = p[3] 
    bb = 1./(bp-1)
    vv = v0/v    
    ee = e0 + b0/bp*v0/vv*(bb*pow(vv,bp) +1.)-bb*b0*v0
    return ee    

def residuals(p,e,v):
    return e - bmeos(v,p)

def bmfit(x,y):
    p = [0,0,0,0] ; p = inifit(x,y)
    return leastsq(residuals, p, args=(y, x))

#-------------------------------------------------------------------------------

def sortstrain(s,e):
    ss=[] ; ee=[] ; ww=[]
    for i in range(len(s)): ww.append(s[i])
    ww.sort()
    for i in range(len(s)):
        ss.append(s[s.index(ww[i])])
        ee.append(e[s.index(ww[i])])
    return ss, ee

#-------------------------------------------------------------------------------

factor     = 2.
startorder = 0

lrydberg = os.path.exists('quantum-espresso') or os.path.exists('vasp')
lplanar  = os.path.exists('planar')
lstartor = os.path.exists('startorder')
lexsinfo = os.path.exists('INFO-elastic-constants')
lexsev   = os.path.exists('energy-vs-volume')

if (lrydberg): factor=1.0

if (lstartor): 
    input_startorder = open('startorder',"r")
    startorder=int(input_startorder.readline().strip().split()[0])

if (not(lexsev)): sys.exit("ERROR: file energy-vs-volume not found!\n")

#-------------------------------------------------------------------------------

bohr_radius     = 0.529177
joule2hartree   = 4.3597482
joule2rydberg   = joule2hartree/2.
unitconv        = joule2hartree/bohr_radius**3*10.**3*factor

electron_charge = 1.602176565e-19 
bohr_radius     = 5.2917721092e-11
rydberg2ev      = 13.605698066      
unitconv        = electron_charge*rydberg2ev/(1e9*bohr_radius**3)*factor

#-------------------------------------------------------------------------------

input_energy = open('energy-vs-volume',"r")
output_birch = open('murnaghan',"w")

energy = [] ; volume = []

while True:
    line = input_energy.readline().strip()
    if len(line) == 0: break
    energy.append(float(line.split()[1]))
    volume.append(float(line.split()[0]))

nv = len(volume)
if (nv < 4): sys.exit("\nERROR: Too few volumes ("+str(nv)+")!\n")

#-------------------------------------------------------------------------------

volume, energy = sortstrain(volume,energy)

#-------------------------------------------------------------------------------

p = bmfit(volume,energy)

v0 = p[0][0] ; b0 = p[0][2]*unitconv ; bp =  p[0][3]

a0sc=(1*p[0][0])**(0.33333333333)
abcc=(2*p[0][0])**(0.33333333333)
afcc=(4*p[0][0])**(0.33333333333)

chi=0 
for i in range(len(volume)): 
    chi=chi+residuals(p[0],energy[i],volume[i])**2
chi=sqrt(chi)/len(volume)

fmt='%10.4f'
amt='%10.4f'
bmt='%8.3f'
pmt='%16.10f'
lmt='%10.2f'

print
print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
print "     V0        B0        BP         a-sc       a-bcc      a-fcc     log(chi)"
print fmt%(v0), bmt%(b0), bmt%(bp)," ",  
print amt%(a0sc), amt%(abcc), amt%(afcc), lmt%(log10(chi))
print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
print
print >> output_birch, pmt%(v0), pmt%(b0), pmt%(bp),  
print >> output_birch, pmt%(a0sc), pmt%(abcc), pmt%(afcc),
print >> output_birch, pmt%(log10(chi))

input_energy.close()
output_birch.close()

#-------------------------------------------------------------------------------








