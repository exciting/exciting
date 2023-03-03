#!/usr/bin/env python2 
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!%!% ----------------------------- ElaStic_Result_Energy_3rd ----------------------------- %!%!%#
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
# python ElaStic_Result_Energy_3rd.py
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
ToGPa = (_e*Ryd2eV)/(1e9*Bohr**3)  # Ryd/[Bohr^3] to GPa
#--------------------------------------------------------------------------------------------------

#%!%!%--- DICTIONARIES ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
head = {\
'CI':'\
    for, space-group number between 207 and 230, Cubic I structure.               \n\n\
    %%%%%%%% 1-ij %%%%%%%%      %%%%%%%% 2-ij %%%%%%%%      %%%%%%%% 3-ij %%%%%%%%  \n\
    111 112 112  0   0   0      112 112 123  0   0   0      112 123 112  0   0   0  \n\
    112 112 123  0   0   0      112 111 112  0   0   0      123 112 112  0   0   0  \n\
    112 123 112  0   0   0      123 112 112  0   0   0      112 112 111  0   0   0  \n\
     0   0   0  144  0   0       0   0   0  155  0   0       0   0   0  155  0   0  \n\
     0   0   0   0  155  0       0   0   0   0  144  0       0   0   0   0  155  0  \n\
     0   0   0   0   0  155      0   0   0   0   0  155      0   0   0   0   0  144 \n\
                                                                                    \n\
    %%%%%%%% 4-ij %%%%%%%%      %%%%%%%% 5-ij %%%%%%%%      %%%%%%%% 6-ij %%%%%%%%  \n\
     0   0   0  144  0   0       0   0   0   0  155  0       0   0   0   0   0  155 \n\
     0   0   0  155  0   0       0   0   0   0  144  0       0   0   0   0   0  155 \n\
     0   0   0  155  0   0       0   0   0   0  155  0       0   0   0   0   0  144 \n\
    144 155 155  0   0   0       0   0   0   0   0  456      0   0   0   0  456  0  \n\
     0   0   0   0   0  456     155 144 155  0   0   0       0   0   0  456  0   0  \n\
     0   0   0   0  456  0       0   0   0  456  0   0      155 155 144  0   0   0  \n',\
'CII':'\
    for, space-group number between 195 and 206, Cubic II structure.              \n\n\
    %%%%%%%% 1-ij %%%%%%%%      %%%%%%%% 2-ij %%%%%%%%      %%%%%%%% 3-ij %%%%%%%%  \n\
    111 112 113  0   0   0      112 112 123  0   0   0      112 123 112  0   0   0  \n\
    112 113 123  0   0   0      112 111 112  0   0   0      123 112 112  0   0   0  \n\
    113 123 112  0   0   0      123 112 113  0   0   0      112 112 111  0   0   0  \n\
     0   0   0  144  0   0       0   0   0  166  0   0       0   0   0  155  0   0  \n\
     0   0   0   0  155  0       0   0   0   0  144  0       0   0   0   0  166  0  \n\
     0   0   0   0   0  166      0   0   0   0   0  155      0   0   0   0   0  144 \n\
                                                                                    \n\
    %%%%%%%% 4-ij %%%%%%%%      %%%%%%%% 5-ij %%%%%%%%      %%%%%%%% 6-ij %%%%%%%%  \n\
     0   0   0  144  0   0       0   0   0   0  155  0       0   0   0   0   0  155 \n\
     0   0   0  155  0   0       0   0   0   0  144  0       0   0   0   0   0  155 \n\
     0   0   0  155  0   0       0   0   0   0  155  0       0   0   0   0   0  144 \n\
    144 155 155  0   0   0       0   0   0   0   0  456      0   0   0   0  456  0  \n\
     0   0   0   0   0  456     155 144 155  0   0   0       0   0   0  456  0   0  \n\
     0   0   0   0  456  0       0   0   0  456  0   0      155 155 144  0   0   0  \n',\
'HI':'\
    for, space-group number between 177 and 194, Hexagonal I structure.           \n\n\
    %%%%%%%% 1-ij %%%%%%%%      %%%%%%%% 2-ij %%%%%%%%      %%%%%%%% 3-ij %%%%%%%%  \n\
    111 112 113  0   0   0      112  A  123  0   0   0      113 123 133  0   0   0  \n\
    112  A  123  0   0   0       A  222 113  0   0   0      123 113 133  0   0   0  \n\
    113 123 133  0   0   0      123 113 133  0   0   0      133 133 333  0   0   0  \n\
     0   0   0  144  0   0       0   0   0  155  0   0       0   0   0  344  0   0  \n\
     0   0   0   0  155  0       0   0   0   0  144  0       0   0   0   0  344  0  \n\
     0   0   0   0   0   B       0   0   0   0   0   C       0   0   0   0   0   D  \n\
                                                                                    \n\
    %%%%%%%% 4-ij %%%%%%%%      %%%%%%%% 5-ij %%%%%%%%      %%%%%%%% 6-ij %%%%%%%%  \n\
     0   0   0  144  0   0       0   0   0   0  155  0       0   0   0   0   0   B  \n\
     0   0   0  155  0   0       0   0   0   0  144  0       0   0   0   0   0   C  \n\
     0   0   0  344  0   0       0   0   0   0  344  0       0   0   0   0   0   D  \n\
    144 155 344  0   0   0       0   0   0   0   0   E       0   0   0   0   E   0  \n\
     0   0   0   0   0   E      155 144 344  0   0   0       0   0   0   E   0   0  \n\
     0   0   0   0   E   0       0   0   0   E   0   0       B   C   D   0   0   0  \n\
                                                                                    \n\
    A = 111 - 222 + 112                                                             \n\
    B = 3/4 * 222 - 1/2 * 111 - 1/4 * 112                                           \n\
    C = 1/2 * 111 - 1/4 * 222 - 1/4 * 112                                           \n\
    D = 1/2 * ( 113 - 123 )                                                         \n\
    E = 1/2 * ( 155 - 144 )                                                         \n',\
'HII':'\
    for, space-group number between 168 and 176, Hexagonal II structure.          \n\n\
    %%%%%%%% 1-ij %%%%%%%%      %%%%%%%% 2-ij %%%%%%%%      %%%%%%%% 3-ij %%%%%%%%  \n\
    111 112 113  0   0  116     112  A  123  0   0   0      113 123 133  0   0   0  \n\
    112  A  123  0   0 -116      A  222 113  0   0  116     123 113 133  0   0   0  \n\
    113 123 133  0   0   0      123 113 133  0   0   0      133 133 333  0   0   0  \n\
     0   0   0  144 145  0       0   0   0  155-145  0       0   0   0  344  0   0  \n\
     0   0   0  145 155  0       0   0   0 -145 144  0       0   0   0   0  344  0  \n\
    116-116  0   0   0   B       0  116  0   0   0   C       0   0   0   0   0   D  \n\
                                                                                    \n\
    %%%%%%%%% 4-ij %%%%%%%%     %%%%%%%%% 5-ij %%%%%%%%     %%%%%%%%% 6-ij %%%%%%%% \n\
     0   0   0  144  0   0       0   0   0   0  155  0       0   0   0   0   0   B  \n\
     0   0   0  155  0   0       0   0   0   0  144  0       0   0   0   0   0   C  \n\
     0   0   0  344  0   0       0   0   0   0  344  0       0   0   0   0   0   D  \n\
    144 155 344  0   0  145      0   0   0   0   0   E       0   0   0   0   E   0  \n\
     0   0   0   0   0   E      155 144 344  0   0 -145      0   0   0   E   0   0  \n\
     0   0   0  145  E   0       0   0   0   E -145  0       B   C   D   0   0 -116 \n\
                                                                                    \n\
    A = 111 - 222 + 112                                                             \n\
    B = 3/4 * 222 - 1/2 * 111 - 1/4 * 112                                           \n\
    C = 1/2 * 111 - 1/4 * 222 - 1/4 * 112                                           \n\
    D = 1/2 * ( 113 - 123 )                                                         \n\
    E = 1/2 * ( 155 - 144 )                                                         \n',\
'RI':'\
    for, space-group number between 149 and 167, Rhombohedral I structure.        \n\n\
    %%%%%%%% 1-ij %%%%%%%%      %%%%%%%% 2-ij %%%%%%%%      %%%%%%%% 3-ij %%%%%%%%  \n\
    111 112 113 114  0   0      112  A  123 124  0   0      113 123 133 134  0   0  \n\
    112  A  123 124  0   0       A  222 113  G   0   0      123 113 133 -134 0   0  \n\
    113 123 133 134  0   0      123 113 133 -134 0   0      133 133 333  0   0   0  \n\
    114 124 134 144  0   0      124  G -134 155  0   0      134 -134 0  344  0   0  \n\
     0   0   0   0  155  F       0   0   0   0  144  H       0   0   0   0  344 134 \n\
     0   0   0   0   F   B       0   0   0   0   H   C       0   0   0   0  134  D  \n\
                                                                                    \n\
    %%%%%%%% 4-ij %%%%%%%%      %%%%%%%% 5-ij %%%%%%%%      %%%%%%%% 6-ij %%%%%%%%  \n\
    114 124 134 144  0   0       0   0   0   0  155  F       0   0   0   0   F   B  \n\
    124  G -134 155  0   0       0   0   0   0  144  H       0   0   0   0   H   C  \n\
    134 -134 0  344  0   0       0   0   0   0  344 134      0   0   0   0  134  D  \n\
    144 155 344 444  0   0       0   0   0   0   0   E       0   0   0   0   E  124 \n\
     0   0   0   0 -444  E      155 144 344  0   0   0       F   H  134  E   0   0  \n\
     0   0   0   0   E  124      F   H  134  E   0   0       B   C   D  124  0   0  \n\
                                                                                    \n\
    A = 111 - 222 + 112                                                             \n\
    B = 3/4 * 222 - 1/2 * 111 - 1/4 * 112                                           \n\
    C = 1/2 * 111 - 1/4 * 222 - 1/4 * 112                                           \n\
    D = 1/2 * ( 113 - 123 )                                                         \n\
    E = 1/2 * ( 155 - 144 )                                                         \n\
    F = 1/2 * ( 114 + 3 * 124 )                                                     \n\
    G =-114 - 2 * 124                                                               \n\
    H = 1/2 * ( 114 - 124 )                                                         \n',\
'RII':'\
    for, space-group number between 143 and 148, Rhombohedral II structure.       \n\n\
    %%%%%%%% 1-ij %%%%%%%%      %%%%%%%% 2-ij %%%%%%%%      %%%%%%%% 3-ij %%%%%%%%  \n\
                                                                                    \n',\
'TI':'\
    for, space-group number between 89 and 142, Tetragonal I structure.           \n\n\
    %%%%%%%% 1-ij %%%%%%%%      %%%%%%%% 2-ij %%%%%%%%      %%%%%%%% 3-ij %%%%%%%%  \n\
                                                                                    \n',\
'TII':'\
    for, space-group number between 75 and 88, Tetragonal II structure.           \n\n\
    %%%%%%%% 1-ij %%%%%%%%      %%%%%%%% 2-ij %%%%%%%%      %%%%%%%% 3-ij %%%%%%%%  \n\
                                                                                    \n',\
'O':'\
    for, space-group number between 16 and 74, Orthorhombic structure.            \n\n\
    %%%%%%%% 1-ij %%%%%%%%      %%%%%%%% 2-ij %%%%%%%%      %%%%%%%% 3-ij %%%%%%%%  \n\
                                                                                    \n',\
'M':'\
    for, space-group number between 3 and 15, Monoclinic structure.               \n\n\
    %%%%%%%% 1-ij %%%%%%%%      %%%%%%%% 2-ij %%%%%%%%      %%%%%%%% 3-ij %%%%%%%%  \n\
                                                                                    \n',\
'N':'\
    for, 1 and 2 space-group number, Triclinic structure.                         \n\n\
    %%%%%%%% 1-ij %%%%%%%%      %%%%%%%% 2-ij %%%%%%%%      %%%%%%%% 3-ij %%%%%%%%  \n\
                                                                                    \n'}
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
    ECs= 56  

elif(3 <= SGN and SGN <= 15):    # Monoclinic
    LC = 'M'
    ECs= 32 

elif(16 <= SGN and SGN <= 74):   # Orthorhombic
    LC = 'O'
    ECs= 20 

elif(75 <= SGN and SGN <= 88):   # Tetragonal II
    LC = 'TII'
    ECs= 16
  
elif(89 <= SGN and SGN <= 142):  # Tetragonal I
    LC = 'TI'
    ECs= 12  

elif(143 <= SGN and SGN <= 148): # Rhombohedral II 
    LC = 'RII'
    ECs= 20

elif(149 <= SGN and SGN <= 167): # Rhombohedral I
    LC = 'RI'
    ECs= 14

elif(168 <= SGN and SGN <= 176): # Hexagonal II
    LC = 'HII'
    ECs= 12

elif(177 <= SGN and SGN <= 194): # Hexagonal I
    LC = 'HI'
    ECs= 10

elif(195 <= SGN and SGN <= 206): # Cubic II
    LC = 'CII'
    ECs=  8

elif(207 <= SGN and SGN <= 230): # Cubic I
    LC = 'CI'
    ECs=  6

else: sys.exit('\n.... Oops ERROR: WRONG Space-Group Number !?!?!?    \n')
#--------------------------------------------------------------------------------------------------

lineseparator=' '
for i in range(0,79):
    lineseparator=lineseparator+'%'

#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!%!% ------------ Calculating the second derivative and Cross-Validation Error ----------- %!%!%#
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
OBJ = open('ElaStic_3rd.in', 'r')

# Range of Deformation
RoD = OBJ.read().strip().split()

if (len(RoD) != 3*ECs):
    sys.exit('\n.... Oops ERROR: Something is WRONG in the "ElaStic_3rd.in" file !?!?!?\n')

os.chdir('Energy-vs-Strain')

A3 = []
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
        sys.exit('\n.... Oops ERROR: NOT enough energy points in "'+Dstn+'-Energy.dat"'\
                 '\n                 for '+str(ordri)+' order polynomial fit.\n')

    coeffs = np.polyfit(strain, energy, ordri)
    A3.append(coeffs[ordri-3])

A3 = np.array(A3)
if (len(A3) != ECs):
    sys.exit('\n.... Oops ERROR: The number of data in the "ElaStic_3rd.in" is NOT equal to '+\
    str(ECs)+'\n')

C = zeros((6,6))

if (cod == 'WIEN2k'):
    CONV = ToGPa * 1.
if (cod == 'exciting'):
    CONV = ToGPa * 2.
if (cod == 'ESPRESSO'):
    CONV = ToGPa * 1.

A3 = A3 * CONV/V0
#%!%!%--- Cubic I structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (LC == 'CI'):
    C111 =-3.*(A3[0] - 4.*A3[1] + A3[3])/2.
    C112 = 1.*(A3[0] - 2.*A3[1] + A3[3])/2.
    C123 = 1.*(A3[0] - 3.*A3[3])/4.
    C144 = 1.*(A3[0] - 4.*A3[1] + A3[3] + 4.*A3[4])/8. 
    C155 = 1.*(A3[0] - 4.*A3[1] + A3[3] + 4.*A3[5])/8. 
    C456 = 1.* A3[2]/8.
#--------------------------------------------------------------------------------------------------

#%!%!%--- Cubic II structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
if (LC == 'CII'):
    C111 =-3.*(A3[0] - 4.*A3[1] + 1.*A3[3])/2.
    C112 = 1.*(A3[0] - 2.*A3[3] + 1.*A3[7])/2.
    C113 = 1.*(A3[0] -12.*A3[1] +10.*A3[3] +1.*A3[7])/2.
    C123 = 1.*(A3[0] + 4.*A3[1] - 4.*A3[3] -3.*A3[7])/4. 
    C144 = 1.*(A3[0] - 4.*A3[1] + 1.*A3[3] +4.*A3[5])/8.
    C155 = 1.*(A3[0] - 4.*A3[1] + 1.*A3[3] +4.*A3[5])/8. 
    C166 = 1.*(A3[0] - 4.*A3[1] + 1.*A3[3] -8.*A3[5] + 12.*A3[6])/8.
    C456 = 1.* A3[2]/8.
#--------------------------------------------------------------------------------------------------

#%!%!%--- Hexagonal I structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (LC == 'HI'):
    C111 = 6.*A3[7]
    C112 = 2.*A3[3] + A3[6] - 4.*A3[7]
    C113 = 1./3.*A3[0] - 1./6.*A3[1] - A3[2] + A3[3] - 0.5*A3[5] - 0.5*A3[6] + A3[7]
    C123 = 1./3.*A3[0] - 7./6.*A3[1] - A3[2] - A3[3] + 0.5*A3[5] - A3[7]
    C133 = 1./3.*A3[0] + 4./3.*A3[1] + A3[2] - 0.5*A3[6]
    C144 =-0.5*(A3[7] - A3[8])
    C155 =-0.5*(A3[3] - A3[9])
    C222 = 6.*A3[3]
    C333 = 6.*A3[2] 
    C344 =-0.5*(A3[2] - A3[4])
#--------------------------------------------------------------------------------------------------

#%!%!%--- Hexagonal II structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
if (LC == 'HII'):
    C111 = 6.*A3[7]
    C112 = 2.*A3[3] + A3[6] - 4.*A3[7]
    C113 = 1./3.*A3[0] - 1./6.*A3[1] - A3[2] - 6.*A3[3] - 0.5*A3[5] - 0.0625*A3[6] + 3.625*A3[7] + 0.875*A3[10]
    C116 = 12.*A3[3] - 0.75*A3[6] - 4.5*A3[7] - 1.5*A3[10]
    C123 = 1./3.*A3[0] - 7./6.*A3[1] - A3[2] + 6.*A3[3] + 0.5*A3[5] - 0.4375*A3[6] - 3.625*A3[7] - 0.875*A3[10]
    C133 = 1./3.*A3[0] + 4./3.*A3[1] + A3[2] - 0.5*A3[6]
    C144 =-0.5*(A3[7] - A3[8])
    C145 = 0.75*(A3[3] - A3[8]- A3[9] + A3[11])
    C155 =-0.50*(A3[3] - A3[9])
    C222 = 6.*A3[3]
    C333 = 6.*A3[2]
    C344 =-0.5*(A3[2] - A3[4])
#--------------------------------------------------------------------------------------------------

#%!%!%--- Rhombohedral I structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
if (LC == 'RI'):
    C111 = 6.*A3[3]
    C112 = A3[1] - 4.*A3[3] + 2.*A3[9]
    C113 =-1./3.*A3[0] + 0.5*A3[1] - 3.*A3[2]- 2.*A3[3] - 4./3.*A3[7] + 2.*A3[8]
    C114 =-2./3.*A3[0] + A3[1] + 2.5*A3[4] - 3.*A3[5] - A3[6] - 2./3.*A3[7] - 4.*A3[9] + 2./3.*A3[11] + 2./3.*A3[12]-0.25*A3[13]
    C123 = A3[0] - A3[1] + A3[2]+ 2.*A3[3] - 2.*A3[8]
    C124 = 1./6.*A3[0] - 0.25*A3[1] - 0.75*A3[4] + A3[5] + 1./6.*A3[7] + A3[9] - 1./6.*A3[11] - 1./6.*A3[12] \
          +0.125*A3[13]
    C133 = 1./3.*A3[0] - 0.5*A3[1] + A3[2]+ 4./3.*A3[7]
    C134 =25./300.*A3[0] - 0.125*A3[1] + 1.5*A3[2] + A3[3] + A3[4] - A3[5] + 1./3.*A3[7] - A3[8] - A3[9] - 0.5*A3[10]+0.25*A3[12]
    C144 = 1./3.*A3[0] - 0.5*A3[1] - 0.5*A3[3] - 1.75*A3[4] + 2.*A3[5] + 0.5*A3[6] + 1./3.*A3[7] + 2.*A3[9] - 1./3.*A3[11] \
          -1./3.*A3[12] + 0.125*A3[13]
    C155 =-0.5*A3[3] + 0.5*A3[6]
    C222 = 6.*A3[9]
    C333 = 6.*A3[2]
    C344 = 0.5*A3[10] - 0.5*A3[2] - 0.5*A3[4]
    C444 = 0.75*A3[4]
#--------------------------------------------------------------------------------------------------

#%!%--- Rhombohedral II structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
#if (LC == 'RII'):
#--------------------------------------------------------------------------------------------------

#%%%--- Tetragonal I structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
#if (LC == 'TI'):
#--------------------------------------------------------------------------------------------------

#%%%--- Tetragonal II structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
#if (LC == 'TII'):
#--------------------------------------------------------------------------------------------------

#%%%--- Orthorhombic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
#if (LC == 'O'):
#--------------------------------------------------------------------------------------------------

#%%%--- Monoclinic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
#if (LC == 'M'):
#--------------------------------------------------------------------------------------------------

#%%%--- Triclinic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
#if (LC == 'N'):
#--------------------------------------------------------------------------------------------------

#%%%--- Write output file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
fo = open('ElaStic_3rd.out','w')
print >>fo,  '    The output of ElaStic code                                            \n'\
             '    Today is '+ time.asctime() +                                         '\n'\
                                                                                       '\n'\
             '    Symmetry of the third-order elastic constants in Voigt notation       \n'\
                + head[LC] +                                                           '\n'

#-- Cubic I Structures ----------------------------------------------------------------------------
if (LC == 'CI'):
    print >>fo, '    C111 =', '%11.2f'%(C111), '  GPa\n'\
                '    C112 =', '%11.2f'%(C112), '  GPa\n'\
                '    C123 =', '%11.2f'%(C123), '  GPa\n'\
                '    C144 =', '%11.2f'%(C144), '  GPa\n'\
                '    C155 =', '%11.2f'%(C155), '  GPa\n'\
                '    C456 =', '%11.2f'%(C456), '  GPa\n'
#--------------------------------------------------------------------------------------------------

#-- Cubic II structures ---------------------------------------------------------------------------
if (LC == 'CII'):
    print >>fo, '    C111 =', '%11.2f'%(C111), '  GPa\n'\
                '    C112 =', '%11.2f'%(C112), '  GPa\n'\
                '    C113 =', '%11.2f'%(C113), '  GPa\n'\
                '    C123 =', '%11.2f'%(C123), '  GPa\n'\
                '    C144 =', '%11.2f'%(C144), '  GPa\n'\
                '    C155 =', '%11.2f'%(C155), '  GPa\n'\
                '    C166 =', '%11.2f'%(C166), '  GPa\n'\
                '    C456 =', '%11.2f'%(C456), '  GPa\n'
#--------------------------------------------------------------------------------------------------

#-- Hexagonal I structures ------------------------------------------------------------------------
if (LC == 'HI'):
    print >>fo, '    C111 =', '%11.2f'%(C111), '  GPa\n'\
                '    C112 =', '%11.2f'%(C112), '  GPa\n'\
                '    C113 =', '%11.2f'%(C113), '  GPa\n'\
                '    C123 =', '%11.2f'%(C123), '  GPa\n'\
                '    C133 =', '%11.2f'%(C133), '  GPa\n'\
                '    C144 =', '%11.2f'%(C144), '  GPa\n'\
                '    C155 =', '%11.2f'%(C155), '  GPa\n'\
                '    C222 =', '%11.2f'%(C222), '  GPa\n'\
                '    C333 =', '%11.2f'%(C333), '  GPa\n'\
                '    C344 =', '%11.2f'%(C344), '  GPa\n'
#--------------------------------------------------------------------------------------------------

#-- Hexagonal II structures -----------------------------------------------------------------------
if (LC == 'HII'):
    print >>fo, '    C111 =', '%11.2f'%(C111), '  GPa\n'\
                '    C112 =', '%11.2f'%(C112), '  GPa\n'\
                '    C113 =', '%11.2f'%(C113), '  GPa\n'\
                '    C116 =', '%11.2f'%(C116), '  GPa\n'\
                '    C123 =', '%11.2f'%(C123), '  GPa\n'\
                '    C133 =', '%11.2f'%(C133), '  GPa\n'\
                '    C144 =', '%11.2f'%(C144), '  GPa\n'\
                '    C145 =', '%11.2f'%(C145), '  GPa\n'\
                '    C155 =', '%11.2f'%(C155), '  GPa\n'\
                '    C222 =', '%11.2f'%(C222), '  GPa\n'\
                '    C333 =', '%11.2f'%(C333), '  GPa\n'\
                '    C344 =', '%11.2f'%(C344), '  GPa\n'
#--------------------------------------------------------------------------------------------------

#-- Rhombohedral I structures ---------------------------------------------------------------------
if (LC == 'RI'):
    print >>fo, '    C111 =', '%11.2f'%(C111), '  GPa\n'\
                '    C112 =', '%11.2f'%(C112), '  GPa\n'\
                '    C113 =', '%11.2f'%(C113), '  GPa\n'\
                '    C114 =', '%11.2f'%(C114), '  GPa\n'\
                '    C123 =', '%11.2f'%(C123), '  GPa\n'\
                '    C124 =', '%11.2f'%(C124), '  GPa\n'\
                '    C133 =', '%11.2f'%(C133), '  GPa\n'\
                '    C134 =', '%11.2f'%(C134), '  GPa\n'\
                '    C144 =', '%11.2f'%(C144), '  GPa\n'\
                '    C155 =', '%11.2f'%(C155), '  GPa\n'\
                '    C222 =', '%11.2f'%(C222), '  GPa\n'\
                '    C333 =', '%11.2f'%(C333), '  GPa\n'\
                '    C344 =', '%11.2f'%(C344), '  GPa\n'\
                '    C444 =', '%11.2f'%(C444), '  GPa\n'
#--------------------------------------------------------------------------------------------------

#--- Rhombohedral II structures -------------------------------------------------------------------
#if (LC == 'RII'):
#--------------------------------------------------------------------------------------------------

#--- Tetragonal I structures ----------------------------------------------------------------------
#if (LC == 'TI'):
#--------------------------------------------------------------------------------------------------

#--- Tetragonal II structures ---------------------------------------------------------------------
#if (LC == 'TII'):
#--------------------------------------------------------------------------------------------------

#--- Orthorhombic structures ----------------------------------------------------------------------
#if (LC == 'O'):
#--------------------------------------------------------------------------------------------------

#--- Monoclinic structures ------------------------------------------------------------------------
#if (LC == 'M'):
#--------------------------------------------------------------------------------------------------

#--- Triclinic structures -------------------------------------------------------------------------
#if (LC == 'N'):
#--------------------------------------------------------------------------------------------------
print >>fo,'\n    ... Have a G00D Day, Week, Month, Year, and Century (if you are lucky) ...    '\
           '\n               Bye-Bye! Tschuess! Ciao! Poka! Zia Jian! KhodaHafez!             \n'
fo.close()
#--------------------------------------------------------------------------------------------------
os.chdir('../')
os.system('cp -f Energy-vs-Strain/ElaStic_3rd.out .')
