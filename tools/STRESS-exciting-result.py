#!/usr/bin/env python 
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!%!% ------------------------------- STRESS-exciting-result ------------------------------ %!%!%#
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#
# AUTHOR:
# Rostam Golesorkhtabar
# r.golesorkhtabar@gmail.com
# 
# DATE:
# Wed Jan 01 00:00:00 2014
#
# SYNTAX:
# python STRESS-exciting-result.py
#
# EXPLANATION:
# 
#__________________________________________________________________________________________________

from sys   import stdin
from numpy import *
import numpy as np
import subprocess
import os.path
import shutil
import math
import time
import sys
import os

#%!%!%--- CONSTANTS ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
_e    =  1.602176565e-19         # elementary charge
Bohr  =  5.291772086e-11         # a.u. to meter
Ha2eV = 27.211396132             # Ha to eV
Tokbar= (_e*Ha2eV)/(1e8*Bohr**3) # Ha/[a.u.]^3 to kbar
ToGPa = (_e*Ha2eV)/(1e9*Bohr**3) # Ha/[a.u.]^3 to GPa
#__________________________________________________________________________________________________

#%!%!%--- Reading the "INFO_Stress" file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
INFO=open('INFO_Stress', 'r')

l1  = INFO.readline()
SGN = int(l1.split()[-1])

l2  = INFO.readline()
V0  = float(l2.split()[-2])

l3  = INFO.readline()
mdr = float(l3.split()[-1])

l4  = INFO.readline()
NoP = int(l4.split()[-1])

if (3 <= SGN and SGN <= 15):
    l5  = INFO.readline()
    unique_axis = str(l5.split()[-1])

INFO.close()
#--------------------------------------------------------------------------------------------------

#%!%!%--- Classify the Space-Group Number ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
if (1 <= SGN and SGN <= 2):      # Triclinic
    LC = 'N'
    SCs= 6

elif(3 <= SGN and SGN <= 15):    # Monoclinic
    LC = 'M'
    SCs= 4

elif(16 <= SGN and SGN <= 74):   # Orthorhombic
    LC = 'O'
    SCs=  3

elif(75 <= SGN and SGN <= 88):   # Tetragonal II
    LC = 'TII'
    SCs=  2
  
elif(89 <= SGN and SGN <= 142):  # Tetragonal I
    LC = 'TI'
    SCs=  2
    
elif(143 <= SGN and SGN <= 148): # Rhombohedral II 
    LC = 'RII'
    SCs=  2
    
elif(149 <= SGN and SGN <= 167): # Rhombohedral I
    LC = 'RI'
    SCs=  2
    
elif(168 <= SGN and SGN <= 176): # Hexagonal II
    LC = 'HII'
    SCs=  2
    
elif(177 <= SGN and SGN <= 194): # Hexagonal I
    LC = 'HI'
    SCs=  2
    
elif(195 <= SGN and SGN <= 206): # Cubic II
    LC = 'CII'
    SCs=  1
    
elif(207 <= SGN and SGN <= 230): # Cubic I
    LC = 'CI'
    SCs=  1
    
else: sys.exit('\n ... Oops ERROR: WRONG Space-Group Number !?!?!?\n')
#--------------------------------------------------------------------------------------------------

lineseparator=' '
for i in range(0,79):
    lineseparator=lineseparator+'%'

#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!%!% ------------ Calculating the first derivative and Cross-Validation Error ------------ %!%!%#
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#

OBJ = open('Stress.in', 'r')

# Range of Deformation
RoD = OBJ.read().strip().split()

if (len(RoD) != 3*SCs):
    sys.exit('\n ... Oops ERROR: Something is WRONG in the "Stress.in" file !?!?!?\n')

os.chdir('Energy-vs-Strain')

A1 = []
for i in range(1, SCs+1):
    if (i<10):
        Dstn = 'Dst0'+str(i)
    else:
        Dstn = 'Dst' +str(i)

    for j in range(0, 3*SCs-1, 3):
        if (RoD[j] == Dstn):
            mdri   = abs(float(RoD[j+1]))
            ordri  = int(abs(float(RoD[j+2])))

    ene_file= open(Dstn+'-Energy.dat', 'r')
    eta_ene = ene_file.read().strip().split()

    strain = []
    energy = []

    for k in range(0, len(eta_ene)-1, 2):
        if (-mdri <= float(eta_ene[k]) and float(eta_ene[k]) <= mdri):
            strain.append(float(eta_ene[k+0]))        
            energy.append(float(eta_ene[k+1]))

    if (len(strain) < ordri+1):  
        sys.exit('\n ... Oops ERROR: NOT enough energy points in "'+Dstn+'-Energy.dat"'
                 '\n                 for '+str(ordri)+' order polynomial fit.\n')

    coeffs = np.polyfit(strain, energy, ordri)
    A1.append(coeffs[ordri-1])

A1 = np.array(A1)
if (len(A1) != SCs):
    sys.exit('\n ... Oops ERROR: The number of data in the "Stress.in" is NOT equal to ' + \
    str(SCs)+'\n')

S = zeros((3,3))

#%!%!%--- Cubic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (LC == 'CI' or \
    LC == 'CII'):
    S[0,0] = A1[0]/3.
    S[1,1] = S[0,0]
    S[2,2] = S[0,0]
#--------------------------------------------------------------------------------------------------

#%!%!%--- Hexagonal, Rhombohedral, and Tetragonal structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (LC == 'HI' or \
    LC == 'HII'or \
    LC == 'RI' or \
    LC == 'RII'or \
    LC == 'TI' or \
    LC == 'TII'):
    S[0,0] = (A1[0] - 1.*A1[1])/3.
    S[1,1] = S[0,0]
    S[2,2] = (A1[0] + 2.*A1[1])/3.
#--------------------------------------------------------------------------------------------------

#%!%!%--- Orthorhombic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
if (LC == 'O'):
    S[0,0] = (A1[0] - 2.*A1[1] - 2.*A1[2])/3. 
    S[1,1] = (A1[0] + 2.*A1[1])/3.
    S[2,2] = (A1[0] + 2.*A1[2])/3.
#--------------------------------------------------------------------------------------------------

#%!%!%--- Monoclinic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
if (LC == 'M'):
    S[0,0] = (A1[0] - 2.*A1[1] - 2.*A1[2])/3. 
    S[1,1] = (A1[0] + 2.*A1[1])/3.
    S[2,2] = (A1[0] + 2.*A1[2])/3.

    if (unique_axis == 'a'): S[1,2] = A1[3]/2.
    if (unique_axis == 'b'): S[0,2] = A1[3]/2.
    if (unique_axis == 'c'): S[0,1] = A1[3]/2.
#--------------------------------------------------------------------------------------------------

#%!%!%--- Triclinic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (LC == 'N'):
    S[0,0] = (A1[0] + 2.* A1[2])/3.,
    S[1,1] = (A1[0] + 2.* A1[1])/3.,
    S[2,2] = (A1[0] - 2.*(A1[1] + A1[2]))/3.,
    S[1,2] = A1[3]/2.,
    S[0,2] = A1[4]/2.,
    S[0,1] = A1[5]/2.
#--------------------------------------------------------------------------------------------------

#%!%!%--- Calculating the Pressure ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%

for i in range(2):
    for j in range(i+1, 3):
        S[j,i] = S[i,j] 

S = S / V0
P = (S[0,0]+S[1,1]+S[2,2])/-3.
#--------------------------------------------------------------------------------------------------

#%!%!%--- Writing the output file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
fo = open('STRESS.OUT', 'w')

print >>fo, ' '+ time.asctime() + '\n'\
          '\n The hydrostatic pressure and stress tensor in kilobar (kbar) \n'\
            ' P = ', round(P*Tokbar, 6), '[kbar]'

for i in range(0,3):
    print >>fo, '',
    for j in range(0,3):
        print >>fo, '%18.6f'%(S[i,j]*Tokbar),
    print >>fo

print >>fo, '\n The hydrostatic pressure and stress tensor in gigapascal (GPa) \n'\
            ' P = ', round(P*ToGPa, 6), '[GPa]'

for i in range(0,3):
    print >>fo, '',
    for j in range(0,3):
        print >>fo, '%18.6f'%(S[i,j]*ToGPa),
    print >>fo

fo.close()
#--------------------------------------------------------------------------------------------------
os.chdir('../')
os.system('cp -f Energy-vs-Strain/STRESS.OUT .')
