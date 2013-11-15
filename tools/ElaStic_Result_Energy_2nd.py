#!/usr/bin/env python 
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!%!% ----------------------------- ElaStic_Result_Energy_2nd ----------------------------- %!%!%#
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#
# AUTHORS:
# Rostam Golesorkhtabar and Pasquale Pavone
# r.golesorkhtabar@gmail.com
# 
# DATE:
# Tue Jan 01 00:00:00 2013
#
# SYNTAX:
# python ElaStic_Result_Energy_2nd.py
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
_e    = 1.602176565e-19            # elementary charge
Bohr  = 5.291772086e-11            # Bohr to meter
Ryd2eV= 13.605698066               # Ryd to eV
ToGPa = (_e*Ryd2eV)/(1e9*Bohr**3)  # Ryd/[a.u.^3] to GPa
#--------------------------------------------------------------------------------------------------

#%!%!%--- DICTIONARIES ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
head = {                                                                     \
'CI':'\
    for, space group-number between 207 and 230, Cubic I structure.        \n\n\
               C11     C12     C12      0       0       0                  \n\
               C12     C11     C12      0       0       0                  \n\
               C12     C12     C11      0       0       0                  \n\
                0       0       0      C44      0       0                  \n\
                0       0       0       0      C44      0                  \n\
                0       0       0       0       0      C44                 \n',\
'CII':'\
    for, space group-number between 195 and 206, Cubic II structure.       \n\n\
               C11     C12     C12      0       0       0                  \n\
               C12     C11     C12      0       0       0                  \n\
               C12     C12     C11      0       0       0                  \n\
                0       0       0      C44      0       0                  \n\
                0       0       0       0      C44      0                  \n\
                0       0       0       0       0      C44                 \n',\
'HI':'\
    for, space group-number between 177 and 194, Hexagonal I structure.    \n\n\
               C11     C12     C13      0       0       0                  \n\
               C12     C11     C13      0       0       0                  \n\
               C13     C13     C33      0       0       0                  \n\
                0       0       0      C44      0       0                  \n\
                0       0       0       0      C44      0                  \n\
                0       0       0       0       0   (C11-C12)/2            \n',\
'HII':'\
    for, space group-number between 168 and 176, Hexagonal II structure.   \n\n\
               C11     C12     C13      0       0       0                  \n\
               C12     C11     C13      0       0       0                  \n\
               C13     C13     C33      0       0       0                  \n\
                0       0       0      C44      0       0                  \n\
                0       0       0       0      C44      0                  \n\
                0       0       0       0       0   (C11-C12)/2            \n',\
'RI':'\
    for, space group-number between 149 and 167, Rhombohedral I structure. \n\n\
               C11     C12     C13     C14      0       0                  \n\
               C12     C11     C13    -C14      0       0                  \n\
               C13     C13     C33      0       0       0                  \n\
               C14    -C14      0      C44      0       0                  \n\
                0       0       0       0      C44     C14                 \n\
                0       0       0       0      C14  (C11-C12)/2            \n',\
'RII':'\
    for, space group-number between 143 and 148, Rhombohedral II structure.\n\n\
               C11     C12     C13     C14     C15      0                  \n\
               C12     C11     C13    -C14    -C15      0                  \n\
               C13     C13     C33      0       0       0                  \n\
               C14    -C14      0      C44      0     -C15                 \n\
               C15    -C15      0       0      C44     C14                 \n\
                0       0       0     -C15     C14  (C11-C12)/2            \n',\
'TI':'\
    for, space group-number between 89 and 142, Tetragonal I structure.    \n\n\
               C11     C12     C13      0       0       0                  \n\
               C12     C11     C13      0       0       0                  \n\
               C13     C13     C33      0       0       0                  \n\
                0       0       0      C44      0       0                  \n\
                0       0       0       0      C44      0                  \n\
                0       0       0       0       0      C66                 \n',\
'TII':'\
    for, space group-number between 75 and 88, Tetragonal II structure.    \n\n\
               C11     C12     C13      0       0      C16                 \n\
               C12     C11     C13      0       0     -C16                 \n\
               C13     C13     C33      0       0       0                  \n\
                0       0       0      C44      0       0                  \n\
                0       0       0       0      C44      0                  \n\
               C16    -C16      0       0       0      C66                 \n',\
'O':'\
    for, space group-number between 16 and 74, Orthorhombic structure.     \n\n\
               C11     C12     C13      0       0       0                  \n\
               C12     C22     C23      0       0       0                  \n\
               C13     C23     C33      0       0       0                  \n\
                0       0       0      C44      0       0                  \n\
                0       0       0       0      C55      0                  \n\
                0       0       0       0       0      C66                 \n',\
'M':'\
    for, space group-number between 3 and 15, Monoclinic structure.        \n\n\
               C11     C12     C13      0       0      C16                 \n\
               C12     C22     C23      0       0      C26                 \n\
               C13     C23     C33      0       0      C36                 \n\
                0       0       0      C44     C45      0                  \n\
                0       0       0      C45     C55      0                  \n\
               C16     C26     C36      0       0      C66                 \n',\
'N':'\
    for, space group-number between 1 and 2, Triclinic structure.          \n\n\
               C11     C12     C13     C14      C15    C16                 \n\
               C12     C22     C23     C24      C25    C26                 \n\
               C13     C23     C33     C34      C35    C36                 \n\
               C14     C24     C34     C44      C45    C46                 \n\
               C15     C25     C35     C45      C55    C56                 \n\
               C16     C26     C36     C46      C56    C66                 \n'}
#--------------------------------------------------------------------------------------------------

#%!%!%--- Reading the "INFO_ElaStic" file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
INFO=open('INFO_ElaStic', 'r')

l1  = INFO.readline()
ordr= int(l1.split()[-1])

if (ordr != 2 and ordr != 3):
    sys.exit('\n.... Oops ERROR: The order of the elastic constant is NOT clear !?!?!?'\
             '\n                 Something is WRONG in the "INFO_ElaStic" file.\n')

l2  = INFO.readline()
mthd= l2.split()[-1]

if (mthd != 'Stress' and mthd != 'Energy'):
    sys.exit('\n.... Oops ERROR: The method of the calculation is NOT clear !?!?!?'\
             '\n                 Something is WRONG in the "INFO_ElaStic" file.\n')

l3  = INFO.readline()
cod = l3.split()[-1]

if (cod != 'WIEN2k' and cod != 'exciting' and cod != 'ESPRESSO'):
    sys.exit('\n.... Oops ERROR: The DFT code is NOT clear !?!?!?'\
             '\n                 Something is WRONG in the "INFO_ElaStic" file.\n')

l4  = INFO.readline()
SGN = int(l4.split()[-1])

l5  = INFO.readline()
V0  = float(l5.split()[-2])

l6  = INFO.readline()
mdr = float(l6.split()[-1])

l7  = INFO.readline()
NoP = int(l7.split()[-1])

INFO.close()
#--------------------------------------------------------------------------------------------------

#%!%!%--- Classify the Space-group number ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
if (1 <= SGN and SGN <= 2):      # Triclinic
    LC = 'N'
    ECs= 21

elif(3 <= SGN and SGN <= 15):    # Monoclinic
    LC = 'M'
    ECs= 13

elif(16 <= SGN and SGN <= 74):   # Orthorhombic
    LC = 'O'
    ECs=  9

elif(75 <= SGN and SGN <= 88):   # Tetragonal II
    LC = 'TII'
    ECs=  7
  
elif(89 <= SGN and SGN <= 142):  # Tetragonal I
    LC = 'TI'
    ECs=  6
    
elif(143 <= SGN and SGN <= 148): # Rhombohedral II 
    LC = 'RII'
    ECs=  7
    
elif(149 <= SGN and SGN <= 167): # Rhombohedral I
    LC = 'RI'
    ECs=  6
    
elif(168 <= SGN and SGN <= 176): # Hexagonal II
    LC = 'HII'
    ECs=  5
    
elif(177 <= SGN and SGN <= 194): # Hexagonal I
    LC = 'HI'
    ECs=  5
    
elif(195 <= SGN and SGN <= 206): # Cubic II
    LC = 'CII'
    ECs=  3
    
elif(207 <= SGN and SGN <= 230): # Cubic I
    LC = 'CI'
    ECs=  3
    
else: sys.exit('\n.... Oops ERROR: WRONG Space-Group Number !?!?!?\n')
#--------------------------------------------------------------------------------------------------

lineseparator=' '
for i in range(0,79):
    lineseparator=lineseparator+'%'

#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!%!% ------------ Calculating the second derivative and Cross-Validation Error ----------- %!%!%#
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
OBJ = open('ElaStic_2nd.in', 'r')

# Range of Deformation
RoD = OBJ.read().strip().split()

if (len(RoD) != 3*ECs):
    sys.exit('\n.... Oops ERROR: Something is WRONG in the "ElaStic_2nd.in" file !?!?!?\n')

os.chdir('Energy-vs-Strain')

A2 = []
for i in range(1, ECs+1):
    if (i<10):
        Dstn = 'Dst0'+str(i)
    else:
        Dstn = 'Dst' +str(i)

    for j in range(0, 3*ECs-1, 3):
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
        sys.exit('\n.... Oops ERROR: NOT enough energy points in "'+Dstn+'-Energy.dat"'
                 '\n                 for '+str(ordri)+' order polynomial fit.\n')

    coeffs = np.polyfit(strain, energy, ordri)
    A2.append(coeffs[ordri-2])

A2 = np.array(A2)
if (len(A2) != ECs):
    sys.exit('\n.... Oops ERROR: The number of data in the "ElaStic_2nd.in" is NOT equal to ' + \
    str(ECs)+'\n')

C = zeros((6,6))

#%!%!%--- Cubic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (LC == 'CI' or \
    LC == 'CII'):
    C[0,0] =-2.*(A2[0]-3.*A2[1])/3.
    C[1,1] = C[0,0]
    C[2,2] = C[0,0]
    C[3,3] = A2[2]/6.
    C[4,4] = C[3,3]
    C[5,5] = C[3,3]
    C[0,1] = (2.*A2[0]-3.*A2[1])/3.
    C[0,2] = C[0,1]
    C[1,2] = C[0,1]
#--------------------------------------------------------------------------------------------------

#%!%!%--- Hexagonal structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (LC == 'HI' or \
    LC == 'HII'):
    C[0,0] = 2.*A2[3]
    C[0,1] = 2./3.*A2[0] + 4./3.*A2[1] - 2.*A2[2] - 2.*A2[3]
    C[0,2] = 1./6.*A2[0] - 2./3.*A2[1] + 0.5*A2[2]
    C[1,1] = C[0,0]
    C[1,2] = C[0,2]
    C[2,2] = 2.*A2[2]
    C[3,3] =-0.5*A2[2] + 0.5*A2[4]
    C[4,4] = C[3,3]
    C[5,5] = .5*(C[0,0] - C[0,1])
#--------------------------------------------------------------------------------------------------

#%!%!%--- Rhombohedral I structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
if (LC == 'RI'):
    C[0,0] = 2.*A2[3]
    C[0,1] = A2[1]- 2.*A2[3]
    C[0,2] = .5*( A2[0] - A2[1] - A2[2])
    C[0,3] = .5*(-A2[3] - A2[4] + A2[5])
    C[1,1] = C[0,0]
    C[1,2] = C[0,2]
    C[1,3] =-C[0,3]
    C[2,2] = 2.*A2[2]
    C[3,3] = .5*A2[4]
    C[4,4] = C[3,3]
    C[4,5] = C[0,3]
    C[5,5] = .5*(C[0,0] - C[0,1])
#--------------------------------------------------------------------------------------------------

#%!%!%--- Rhombohedral II structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (LC == 'RII'):
    C[0,0] = 2.*A2[3]
    C[0,1] = A2[1]- 2.*A2[3]
    C[0,2] = .5*( A2[0] - A2[1] - A2[2])
    C[0,3] = .5*(-A2[3] - A2[4] + A2[5])
    C[0,4] = .5*(-A2[3] - A2[4] + A2[6])
    C[1,1] = C[0,0]
    C[1,2] = C[0,2]
    C[1,3] =-C[0,3]
    C[1,4] =-C[0,4]    
    C[2,2] = 2.*A2[2]
    C[3,3] = .5*A2[4]
    C[3,5] =-C[0,4]
    C[4,4] = C[3,3]
    C[4,5] = C[0,3]
    C[5,5] = .5*(C[0,0] - C[0,1])
#--------------------------------------------------------------------------------------------------

#%!%!%--- Tetragonal I structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
if (LC == 'TI'):
    C[0,0] = (A2[0]+2.*A2[1])/3.+.5*A2[2]-A2[3]
    C[0,1] = (A2[0]+2.*A2[1])/3.-.5*A2[2]-A2[3]
    C[0,2] = A2[0]/6.-2.*A2[1]/3.+.5*A2[3]
    C[1,1] = C[0,0]
    C[1,2] = C[0,2]
    C[2,2] = 2.*A2[3]
    C[3,3] = .5*A2[4]
    C[4,4] = C[3,3]
    C[5,5] = .5*A2[5]
#--------------------------------------------------------------------------------------------------

#%!%!%--- Tetragonal II structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (LC == 'TII'):
    C[0,0] = (A2[0]+2.*A2[1])/3.+.5*A2[2]-A2[4]
    C[1,1] = C[0,0]
    C[0,1] = (A2[0]+2.*A2[1])/3.-.5*A2[2]-A2[4]
    C[0,2] = A2[0]/6.-(2./3.)*A2[1]+.5*A2[4]
    C[0,5] = (-A2[2]+A2[3]-A2[6])/4.
    C[1,2] = C[0,2]
    C[1,5] =-C[0,5]
    C[2,2] = 2.*A2[4]
    C[3,3] = .5*A2[5]
    C[4,4] = C[3,3]
    C[5,5] = .5*A2[6]
#--------------------------------------------------------------------------------------------------

#%!%!%--- Orthorhombic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
if (LC == 'O'):
    C[0,0] = 2.*A2[0]/3.+4.*A2[1]/3.+A2[3]-2.*A2[4]-2.*A2[5]
    C[0,1] = 1.*A2[0]/3.+2.*A2[1]/3.-.5*A2[3]-A2[5]
    C[0,2] = 1.*A2[0]/3.-2.*A2[1]/3.+4.*A2[2]/3.-.5*A2[3]-A2[4]
    C[1,1] = 2.*A2[4]
    C[1,2] =-2.*A2[1]/3.-4.*A2[2]/3.+.5*A2[3]+A2[4]+A2[5]
    C[2,2] = 2.*A2[5]
    C[3,3] = .5*A2[6]
    C[4,4] = .5*A2[7]
    C[5,5] = .5*A2[8]
#--------------------------------------------------------------------------------------------------

#%!%!%--- Monoclinic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
if (LC == 'M'):
    C[0,0] = 2.*A2[0]/3.+8.*(A2[1]+A2[2])/3.-2.*(A2[5]+A2[8]+A2[9])
    C[0,1] = A2[0]/3.+4.*(A2[1]+A2[2])/3.-2.*A2[5]-A2[9]
    C[0,2] =(A2[0]-4.*A2[2])/3.+A2[5]-A2[8]
    C[0,5] =-1.*A2[0]/6.-2.*(A2[1]+A2[2])/3.+.5*(A2[5]+A2[7]+A2[8]+A2[9]-A2[12])
    C[1,1] = 2.*A2[8]
    C[1,2] =-4.*(2.*A2[1]+A2[2])/3.+2.*A2[5]+A2[8]+A2[9]+A2[12]
    C[1,5] =-1.*A2[0]/6.-2.*(A2[1]+A2[2])/3.-.5*A2[3]+A2[5]+.5*(A2[7]+A2[8]+A2[9])
    C[2,2] = 2.*A2[9]
    C[2,5] =-1.*A2[0]/6.+2.*A2[1]/3.-.5*(A2[3]+A2[4]-A2[7]-A2[8]-A2[9]-A2[12])
    C[3,3] = .5*A2[10]
    C[3,4] = .25*(A2[6]-A2[10]-A2[11])
    C[4,4] = .5*A2[11]
    C[5,5] = .5*A2[12]
#--------------------------------------------------------------------------------------------------

#%!%!%--- Triclinic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (LC == 'N'):
    C[0,0] = 2.*A2[0]
    C[0,1] = 1.*(-A2[0]-A2[1]+A2[6])
    C[0,2] = 1.*(-A2[0]-A2[2]+A2[7])
    C[0,3] = .5*(-A2[0]-A2[3]+A2[8]) 
    C[0,4] = .5*(-A2[0]+A2[9]-A2[4])
    C[0,5] = .5*(-A2[0]+A2[10]-A2[5])
    C[1,1] = 2.*A2[1]
    C[1,2] = 1.*(A2[11]-A2[1]-A2[2])
    C[1,3] = .5*(A2[12]-A2[1]-A2[3])
    C[1,4] = .5*(A2[13]-A2[1]-A2[4])
    C[1,5] = .5*(A2[14]-A2[1]-A2[5])
    C[2,2] = 2.*A2[2] 
    C[2,3] = .5*(A2[15]-A2[2]-A2[3])
    C[2,4] = .5*(A2[16]-A2[2]-A2[4])
    C[2,5] = .5*(A2[17]-A2[2]-A2[5])
    C[3,3] = .5*A2[3]
    C[3,4] = .25*(A2[18]-A2[3]-A2[4])
    C[3,5] = .25*(A2[19]-A2[3]-A2[5])
    C[4,4] = .5*A2[4]
    C[4,5] = .25*(A2[20]-A2[4]-A2[5])
    C[5,5] = .5*A2[5]
#--------------------------------------------------------------------------------------------------

#%!%!%--- Calculating the elastic moduli ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (cod == 'WIEN2k'):
    CONV = ToGPa * 1.
if (cod == 'exciting'):
    CONV = ToGPa * 2.
if (cod == 'ESPRESSO'):
    CONV = ToGPa * 1.

for i in range(5):
    for j in range(i+1,6):
        C[j,i] = C[i,j] 

C = C * CONV/V0

BV = (C[0,0]+C[1,1]+C[2,2]+2*(C[0,1]+C[0,2]+C[1,2]))/9
GV = ((C[0,0]+C[1,1]+C[2,2])-(C[0,1]+C[0,2]+C[1,2])+3*(C[3,3]+C[4,4]+C[5,5]))/15
EV = (9*BV*GV)/(3*BV+GV)
nuV= (1.5*BV-GV)/(3*BV+GV)
S  = linalg.inv(C)
BR = 1/(S[0,0]+S[1,1]+S[2,2]+2*(S[0,1]+S[0,2]+S[1,2]))
GR =15/(4*(S[0,0]+S[1,1]+S[2,2])-4*(S[0,1]+S[0,2]+S[1,2])+3*(S[3,3]+S[4,4]+S[5,5]))
ER = (9*BR*GR)/(3*BR+GR)
nuR= (1.5*BR-GR)/(3*BR+GR)
BH = 0.50*(BV+BR)
GH = 0.50*(GV+GR)
EH = (9.*BH*GH)/(3.*BH+GH)
nuH= (1.5*BH-GH)/(3.*BH+GH)
AVR= 100.*(GV-GR)/(GV+GR)
#--------------------------------------------------------------------------------------------------

#%!%!%--- Writing the output file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
fo = open('ElaStic_2nd.out', 'w')
print >>fo,  '    The output of ElaStic code                                              \n'\
             '    Today is '+ time.asctime() +                                           '\n'\
                                                                                         '\n'\
             '    Symmetry of the second-order elastic constant matrix in Voigt notation. \n'\
                + head[LC] +                                                             '\n'\
             '    Elastic constant (stiffness) matrix in GPa:                             \n'

for i in range(0,6):
    print >>fo, '',
    for j in range(0,6):
        print >>fo, '%11.1f'%(C[i,j]),
    print >>fo

print >>fo,'\n\n    Elastic compliance matrix in 1/GPa: \n'

for i in range(0,6):
    print >>fo, '',
    for j in range(0,6):
        print >>fo, '%11.5f'%(S[i,j]),
    print >>fo

print >>fo, '\n'+ lineseparator +'\n'

print >>fo, '    Voigt bulk  modulus, B_V = {0}  GPa'.format('%8.2f'%(BV))
print >>fo, '    Voigt shear modulus, G_V = {0}  GPa'.format('%8.2f'%(GV)) + '\n'

print >>fo, '    Reuss bulk  modulus, B_R = {0}  GPa'.format('%8.2f'%(BR))
print >>fo, '    Reuss shear modulus, G_R = {0}  GPa'.format('%8.2f'%(GR)) + '\n'

print >>fo, '    Hill bulk  modulus,  B_H = {0}  GPa'.format('%8.2f'%(BH))
print >>fo, '    Hill shear modulus,  G_H = {0}  GPa'.format('%8.2f'%(GH))

print >>fo, '\n'+ lineseparator +'\n'

print >>fo, '    Voigt Young modulus,  E_V = {0}  GPa'.format('%8.2f'%(EV))
print >>fo, '    Voigt Poisson ratio, nu_V = {0}'     .format('%8.2f'%(nuV)) + '\n'

print >>fo, '    Reuss Young modulus,  E_R = {0}  GPa'.format('%8.2f'%(ER))
print >>fo, '    Reuss Poisson ratio, nu_R = {0}'     .format('%8.2f'%(nuR)) + '\n'

print >>fo, '    Hill Young modulus,   E_H = {0}  GPa'.format('%8.2f'%(EH))
print >>fo, '    Hill Poisson ratio,  nu_H = {0}'     .format('%8.2f'%(nuH))

print >>fo, '\n'+ lineseparator +'\n'

print >>fo, '    Elastic Anisotropy in polycrystalline, AVR = {0} %'.format('%8.3f'%(AVR))

print >>fo, '\n'+ lineseparator +'\n'

print >>fo, '    Eigenvalues of elastic constant (stiffness) matrix:   \n'

eigval=linalg.eig(C)
for i in range(6):
    print >>fo,'%16.1f' % float(eigval[0][i])

print >>fo,'\n    ... Have a G00D Day, Week, Month, Year, and Century (if you are lucky) ...    '\
           '\n               Bye-Bye! Tschuess! Ciao! Poka! Zia Jian! KhodaHafez!             \n'
fo.close()
#--------------------------------------------------------------------------------------------------
os.chdir('../')
os.system('cp -f Energy-vs-Strain/ElaStic_2nd.out .')
