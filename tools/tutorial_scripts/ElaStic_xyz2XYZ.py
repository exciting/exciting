#!/usr/bin/env python 
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!%!% ---------------------------------- ElaStic_xyz2XYZ ---------------------------------- %!%!%#
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
# python ElaStic_xyz2XYZ.py 
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

#%!%!%--- SUBROUTINS AND FUNCTIONS ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
def C9nmop(m,n,o,p):
    S=0.
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    S = S + TM[m,i]*TM[n,j]*TM[o,k]*TM[p,l]*C9o[i,j,k,l]
    C9n[m,n,o,p] = S
    return C9n[m,n,o,p]
#--------------------------------------------------------------------------------------------------

def myabs(x):
    y=sqrt(x**2)
    return y
#--------------------------------------------------------------------------------------------------

#%!%!%--- Checking the Transformation_Matrix file exist ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
if (os.path.exists('Transformation_Matrix') == False):
    sys.exit('\n.... Oops ERROR: There is NO file with "Transformation_Matrix" name.\n\n\
     In order to transform elastic constant matrix from old cartision coordinate\n\
     System to the new one. ElaStic needs two input files:                      \n\
                                                                                \n\
--=> 1- The output file of ElaStic program.                                     \n\
--=> 2- Transformation_Matrix file, this file contains 3x3 matrix like below:   \n\n\
                           i\'.i   i\'.j   i\'.k                                \n\
                           j\'.i   j\'.j   j\'.k                                \n\
                           k\'.i   k\'.j   k\'.k                                \n\n\
     i\', j\' and k\' are unit vectors in new coordinate systems and            \n\
     i, j and k are unit vectors in old coordinate systems.                     \n')

#--------------------------------------------------------------------------------------------------

lineseparator=' '
for i in range(0,79):
    lineseparator=lineseparator+'%'

#%!%!%--- Reading the Transformation_Matrix file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
TM_f = open('Transformation_Matrix', 'r')
TM = zeros((3,3))

for i in range(3):
    line = TM_f.readline()
    for j in range(3):
        TM[i,j] = line.split()[j]

if (myabs(myabs(np.linalg.det(TM))-1.) > 1e-10):
    sys.exit('\n.... Oops ERROR: TRANSFORMATION MATRIX IS WRONG'\
             '\n                 Determinant of matrix should be one.\n')

for i in range(3):
    sumsq = TM[i,0]**2 + TM[i,1]**2 + TM[i,2]**2
    if (myabs(myabs(sumsq)-1.) > 1e-10): 
       sys.exit('\n.... Oops ERROR: TRANSFORMATION MATRIX IS WRONG'\
                '\n                 Sum of squares on each rows should be one.\n')  

#--------------------------------------------------------------------------------------------------
infile = 'ElaStic_2nd.out'
outfile= 'xyz2XYZ_2nd.out'

if (os.path.exists(infile)):
    inf = open(infile, 'r')
    for i in range(15):
        inf.readline()
        
    C6o = zeros((6,6))
    for i in range(6):
        line = inf.readline()
        for j in range(6):
            C6o[i,j] = line.split()[j]

    C9o = zeros((3,3,3,3))
    #%!%!%--- Making the C9o Matrix ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
    C9o[0,0,0,0] = C6o[0,0]; C9o[0,0,0,1] = C6o[0,5]; C9o[0,0,0,2] = C6o[0,4]
    C9o[0,0,1,0] = C6o[0,5]; C9o[0,0,1,1] = C6o[0,1]; C9o[0,0,1,2] = C6o[0,3]
    C9o[0,0,2,0] = C6o[0,4]; C9o[0,0,2,1] = C6o[0,3]; C9o[0,0,2,2] = C6o[0,2]

    C9o[0,1,0,0] = C6o[0,5]; C9o[0,1,0,1] = C6o[5,5]; C9o[0,1,0,2] = C6o[4,5]
    C9o[0,1,1,0] = C6o[5,5]; C9o[0,1,1,1] = C6o[1,5]; C9o[0,1,1,2] = C6o[3,5]
    C9o[0,1,2,0] = C6o[4,5]; C9o[0,1,2,1] = C6o[3,5]; C9o[0,1,2,2] = C6o[2,5]

    C9o[0,2,0,0] = C6o[0,4]; C9o[0,2,0,1] = C6o[4,5]; C9o[0,2,0,2] = C6o[4,4]
    C9o[0,2,1,0] = C6o[4,5]; C9o[0,2,1,1] = C6o[1,4]; C9o[0,2,1,2] = C6o[3,4]
    C9o[0,2,2,0] = C6o[4,4]; C9o[0,2,2,1] = C6o[3,4]; C9o[0,2,2,2] = C6o[2,4]
    #************************************************************************
    C9o[1,0,0,0] = C6o[0,5]; C9o[1,0,0,1] = C6o[5,5]; C9o[1,0,0,2] = C6o[4,5]
    C9o[1,0,1,0] = C6o[5,5]; C9o[1,0,1,1] = C6o[1,5]; C9o[1,0,1,2] = C6o[3,5]
    C9o[1,0,2,0] = C6o[4,5]; C9o[1,0,2,1] = C6o[3,5]; C9o[1,0,2,2] = C6o[2,5]

    C9o[1,1,0,0] = C6o[0,1]; C9o[1,1,0,1] = C6o[1,5]; C9o[1,1,0,2] = C6o[1,4]
    C9o[1,1,1,0] = C6o[1,5]; C9o[1,1,1,1] = C6o[1,1]; C9o[1,1,1,2] = C6o[1,3]
    C9o[1,1,2,0] = C6o[1,4]; C9o[1,1,2,1] = C6o[1,3]; C9o[1,1,2,2] = C6o[1,2]

    C9o[1,2,0,0] = C6o[0,3]; C9o[1,2,0,1] = C6o[3,5]; C9o[1,2,0,2] = C6o[3,4]
    C9o[1,2,1,0] = C6o[3,5]; C9o[1,2,1,1] = C6o[1,3]; C9o[1,2,1,2] = C6o[3,3]
    C9o[1,2,2,0] = C6o[3,4]; C9o[1,2,2,1] = C6o[3,3]; C9o[1,2,2,2] = C6o[2,3]
    #************************************************************************
    C9o[2,0,0,0] = C6o[0,4]; C9o[2,0,0,1] = C6o[4,5]; C9o[2,0,0,2] = C6o[4,4]
    C9o[2,0,1,0] = C6o[4,5]; C9o[2,0,1,1] = C6o[1,4]; C9o[2,0,1,2] = C6o[3,4]
    C9o[2,0,2,0] = C6o[4,4]; C9o[2,0,2,1] = C6o[3,4]; C9o[2,0,2,2] = C6o[2,4]

    C9o[2,1,0,0] = C6o[0,3]; C9o[2,1,0,1] = C6o[3,5]; C9o[2,1,0,2] = C6o[3,4]
    C9o[2,1,1,0] = C6o[3,5]; C9o[2,1,1,1] = C6o[1,3]; C9o[2,1,1,2] = C6o[3,3]
    C9o[2,1,2,0] = C6o[3,4]; C9o[2,1,2,1] = C6o[3,3]; C9o[2,1,2,2] = C6o[2,3]

    C9o[2,2,0,0] = C6o[0,2]; C9o[2,2,0,1] = C6o[2,5]; C9o[2,2,0,2] = C6o[2,4]
    C9o[2,2,1,0] = C6o[2,5]; C9o[2,2,1,1] = C6o[1,2]; C9o[2,2,1,2] = C6o[2,3]
    C9o[2,2,2,0] = C6o[2,4]; C9o[2,2,2,1] = C6o[2,3]; C9o[2,2,2,2] = C6o[2,2]
    #----------------------------------------------------------------------------------------------
    C6n = zeros((6,6))
    C9n = zeros((3,3,3,3))
    C9nmop(0,0,0,0) ;  C6n[0,0] = C9n[0,0,0,0]
    C9nmop(0,0,1,1) ;  C6n[0,1] = C9n[0,0,1,1] ; C6n[1,0] = C6n[0,1] 
    C9nmop(0,0,2,2) ;  C6n[0,2] = C9n[0,0,2,2] ; C6n[2,0] = C6n[0,2]
    C9nmop(0,0,1,2) ;  C6n[0,3] = C9n[0,0,1,2] ; C6n[3,0] = C6n[0,3]
    C9nmop(0,0,0,2) ;  C6n[0,4] = C9n[0,0,0,2] ; C6n[4,0] = C6n[0,4]
    C9nmop(0,0,0,1) ;  C6n[0,5] = C9n[0,0,0,1] ; C6n[5,0] = C6n[0,5]

    C9nmop(1,1,1,1) ;  C6n[1,1] = C9n[1,1,1,1]
    C9nmop(1,1,2,2) ;  C6n[1,2] = C9n[1,1,2,2] ; C6n[2,1] = C6n[1,2]
    C9nmop(1,1,1,2) ;  C6n[1,3] = C9n[1,1,1,2] ; C6n[3,1] = C6n[1,3]
    C9nmop(1,1,0,2) ;  C6n[1,4] = C9n[1,1,0,2] ; C6n[4,1] = C6n[1,4]
    C9nmop(1,1,0,1) ;  C6n[1,5] = C9n[1,1,0,1] ; C6n[5,1] = C6n[1,5]

    C9nmop(2,2,2,2) ;  C6n[2,2] = C9n[2,2,2,2] 
    C9nmop(2,2,1,2) ;  C6n[2,3] = C9n[2,2,1,2] ; C6n[3,2] = C6n[2,3]
    C9nmop(2,2,0,2) ;  C6n[2,4] = C9n[2,2,0,2] ; C6n[4,2] = C6n[2,4]
    C9nmop(2,2,0,1) ;  C6n[2,5] = C9n[2,2,0,1] ; C6n[5,2] = C6n[2,5]

    C9nmop(1,2,1,2) ;  C6n[3,3] = C9n[1,2,1,2]
    C9nmop(1,2,0,2) ;  C6n[3,4] = C9n[1,2,0,2] ; C6n[4,3] = C6n[3,4]
    C9nmop(1,2,0,1) ;  C6n[3,5] = C9n[1,2,0,1] ; C6n[5,3] = C6n[3,5]

    C9nmop(0,2,0,2) ;  C6n[4,4] = C9n[0,2,0,2]
    C9nmop(0,2,0,1) ;  C6n[4,5] = C9n[0,2,0,1] ; C6n[5,4] = C6n[4,5]

    C9nmop(0,1,0,1) ;  C6n[5,5] = C9n[0,1,0,1]

    #%!%!%--- Calculating the elastic moduli ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
    BV = (C6n[0,0]+C6n[1,1]+C6n[2,2]+2*(C6n[0,1]+C6n[0,2]+C6n[1,2]))/9
    GV = ((C6n[0,0]+C6n[1,1]+C6n[2,2])-(C6n[0,1]+C6n[0,2]+C6n[1,2])+3*(C6n[3,3]+C6n[4,4]+C6n[5,5]))/15
    EV = (9*BV*GV)/(3*BV+GV)
    nuV= (1.5*BV-GV)/(3*BV+GV)
    S  = linalg.inv(C6n)
    BR = 1/(S[0,0]+S[1,1]+S[2,2]+2*(S[0,1]+S[0,2]+S[1,2]))
    GR =15/(4*(S[0,0]+S[1,1]+S[2,2])-4*(S[0,1]+S[0,2]+S[1,2])+3*(S[3,3]+S[4,4]+S[5,5]))
    ER = (9*BR*GR)/(3*BR+GR)
    nuR= (1.5*BR-GR)/(3*BR+GR)
    BH = 0.50*(BV+BR)
    GH = 0.50*(GV+GR)
    EH = (9.*BH*GH)/(3.*BH+GH)
    nuH= (1.5*BH-GH)/(3.*BH+GH)
    AVR= 100.*(GV-GR)/(GV+GR) 
    #----------------------------------------------------------------------------------------------

#%!%!%--- Writing the output file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
fo = open(outfile, 'w')
print >>fo,  '    The output of ElaStic_xyz2XYZ code                                    \n'\
             '    Today is '+ time.asctime() +                                         '\n'\
             '                                                                          \n'\
             '                                                                          \n'\
             '    Elastic constant in new cartesian coordinate system                   \n'\
             '    It has been changed based on transformation matrix.                   \n'\
             '                                                                          \n'\
             '    Transformation Matrix:                                                \n'

for i in range(3):
    print >>fo,'    '+str(TM[i,0])+'  '+str(TM[i,1])+'  '+str(TM[i,2])
             
print >>fo, '\n    Elastic constant (stiffness) matrix in GPa:                          \n'

for i in range(0,6):
    print >>fo, '',
    for j in range(0,6):
        print >>fo, '%11.1f'%(C6n[i,j]),
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

eigval=linalg.eig(C6n)
for i in range(6):
    print >>fo,'%16.1f' % float(eigval[0][i])

print >>fo,'\n    ... Have a G00D Day, Week, Month, Year, and Century (if you are lucky) ...    '\
           '\n               Bye-Bye! Tschuess! Ciao! Poka! Zia Jian! KhodaHafez!             \n'
fo.close()
#--------------------------------------------------------------------------------------------------
