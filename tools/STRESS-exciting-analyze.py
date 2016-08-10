#!/usr/bin/env python
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!%!% ------------------------------ STRESS-exciting-analyze ------------------------------ %!%!%#
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
# python STRESS-exciting-analyze.py
# 
# EXPLANATION:
# 
#__________________________________________________________________________________________________

from sys   import stdin
from numpy import *
import numpy as np
import subprocess
import warnings
import os.path
import shutil
import copy
import math
import sys
import os

#%!%!%--- CONSTANTS ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
_e    =  1.602176565e-19         # elementary charge
Bohr  =  5.291772086e-11         # a.u. to meter
Ha2eV = 27.211396132             # Ha to eV
Tokbar= (_e*Ha2eV)/(1e8*Bohr**3) # Ha/[a.u.]^3 to kbar
ToGPa = (_e*Ha2eV)/(1e9*Bohr**3) # Ha/[a.u.]^3 to GPa
#__________________________________________________________________________________________________

#%!%!%--- SUBROUTINS AND FUNCTIONS ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
def sortlist(lst1, lst2):
    temp = copy.copy(lst1)

    lst3 = []
    lst4 = []

    temp.sort()

    for i in range(len(lst1)):
        lst3.append(lst1[lst1.index(temp[i])])
        lst4.append(lst2[lst1.index(temp[i])])

    return lst3, lst4
#--------------------------------------------------------------------------------------------------

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

INFO.close()
#--------------------------------------------------------------------------------------------------

#%!%!%--- Calculating the Space-Group Number and classifying it ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
if (1 <= SGN and SGN <= 2):      # Triclinic
    LC = 'N'
    SCs= 6

elif(3 <= SGN and SGN <= 15):    # Monoclinic
    LC = 'M'
    SCs= 4

elif(16 <= SGN and SGN <= 74):   # Orthorhombic
    LC = 'O'
    SCs= 3

elif(75 <= SGN and SGN <= 88):   # Tetragonal II
    LC = 'TII'
    SCs= 2

elif(89 <= SGN and SGN <= 142):  # Tetragonal I
    LC = 'TI'
    SCs= 2

elif(143 <= SGN and SGN <= 148): # Rhombohedral II 
    LC = 'RII'
    SCs= 2

elif(149 <= SGN and SGN <= 167): # Rhombohedral I
    LC = 'RI'
    SCs= 2

elif(168 <= SGN and SGN <= 176): # Hexagonal II
    LC = 'HII'
    SCs= 2

elif(177 <= SGN and SGN <= 194): # Hexagonal I
    LC = 'HI'
    SCs= 2

elif(195 <= SGN and SGN <= 206): # Cubic II
    LC = 'CII'
    SCs= 1

elif(207 <= SGN and SGN <= 230): # Cubic I
    LC = 'CI'
    SCs= 1

else: sys.exit('\n.... Oops ERROR: WRONG Space-Group Number !?!?!?    \n')
#--------------------------------------------------------------------------------------------------

#%!%!%--- Reading the energies ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
for i in range(1, SCs+1):
    if (i<10):
        Dstn = 'Dst0'+ str(i)
    else: 
        Dstn = 'Dst' + str(i)
    if (os.path.exists(Dstn) == False):
        sys.exit('.... Oops ERROR: Where is the '+ Dstn +' directory !?!?!?    \n')    
    os.chdir(Dstn)

    f = open(Dstn+'-Energy.dat', 'w')

    for j in range(1, NoP+1):
        if (j<10):
            Dstn_num = Dstn +'_0'+str(j)
        else:
            Dstn_num = Dstn +'_' +str(j)

        if (os.path.exists(Dstn_num)):
            os.chdir(Dstn_num)

        if (os.path.exists('INFO.OUT')):
            for line in open('INFO.OUT', 'r'):
                if (line.find('Total energy                               :')>=0): 
                    energy = float(line.split()[-1])

        s = j-(NoP+1)/2
        r = 2*mdr*s/(NoP-1)
        #if (s==0): r=0.0001

        if (r>0):
            strain ='+%12.10f'%r
        else:
            strain = '%13.10f'%r

        if (energy != []):
            print >>f, strain,'   ', energy
        os.chdir('../')

    f.close()
    os.chdir('../')
#--------------------------------------------------------------------------------------------------

warnings.simplefilter('ignore', np.RankWarning)

#%!%!%--- Directory Management ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (os.path.exists('Energy-vs-Strain_old')):
    shutil.rmtree( 'Energy-vs-Strain_old')   

if (os.path.exists('Energy-vs-Strain')):
    os.rename(     'Energy-vs-Strain','Energy-vs-Strain_old')

os.mkdir('Energy-vs-Strain')
os.chdir('Energy-vs-Strain')

os.system('cp -f ../Dst??/Dst??-Energy.dat .')

#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!%!% ------------ Calculating the first derivative and Cross-Validation Error ------------ %!%!%#
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#

for i in range(1, SCs+1):
    if (i<10):
        Dstn = 'Dst0'+str(i)
    else:
        Dstn = 'Dst' +str(i)

    fD = open(Dstn+'_d1E.dat', 'w')

    fE = open(Dstn+'_CVe.dat', 'w')
    print >> fD, '# Max. eta    SUM(sig_i) \n#'
    print >> fE, '# Max. eta    Cross-Validation Error   \n#'

    for j in range(5, 0, -2):
        if (j == 1): nth = '1st'
        if (j == 2): nth = '2nd'
        if (j == 3): nth = '3rd'
        if (j == 4): nth = '4th'
        if (j == 5): nth = '5th'

        print >> fD, '\n# '+ nth +' order fit.'
        print >> fE, '\n# '+ nth +' order fit.'

        #--- Reading the input files --------------------------------------------------------------
        eta_ene= open(Dstn+'-Energy.dat', 'r')

        nl     = 0
        strain = []
        energy = []
        while (nl < NoP):
            line = eta_ene.readline()
            if (line == ''): break 
            line = line.strip().split()
            if (len(line) == 2): 
                nl +=1
                eta, ene = line
                strain.append(float(eta))
                energy.append(float(ene)) 
            elif (len(line) == 0): pass
            else:
                sys.exit('\n ... Oops ERROR: Strain and Energy are NOT defined correctly in "' +\
                          Dstn+'-Energy.dat" !?!?!?\n')

        eta_ene.close()
        strain, energy = sortlist(strain, energy)
        strain0 = copy.copy(strain)
        energy0 = copy.copy(energy)
        #------------------------------------------------------------------------------------------

        #------------------------------------------------------------------------------------------
        while (len(strain) > j): 
            emax  = max(strain)
            emin  = min(strain)
            emax  = max(abs(emin),abs(emax))
            coeffs= polyfit(strain, energy, j)
            sig_i = coeffs[j-1]*Tokbar/V0         # in kbar unit

            print >>fD, '%13.10f'%emax, '%18.6f'%sig_i

            if (abs(strain[0]+emax) < 1.e-7):
                strain.pop(0); energy.pop(0)
            if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                strain.pop()
                energy.pop()

        #--- Cross-Validation error calculations --------------------------------------------------
        strain = copy.copy(strain0)
        energy = copy.copy(energy0)
        while (len(strain) > j+1): 
            emax = max(strain)
            emin = min(strain)
            emax = max(abs(emin),abs(emax))

            S = 0
            for k in range(len(strain)):
                Y      = energy[k]
                etatmp = []
                enetmp = []

                for l in range(len(strain)):
                    if (l==k): pass
                    else:            
                        etatmp.append(strain[l])
                        enetmp.append(energy[l])

                Yfit = polyval(polyfit(etatmp,enetmp, j), strain[k])
                S    = S + (Yfit-Y)**2

            CV = sqrt(S/len(strain))
            print >>fE, '%13.10f'%emax, CV

            if (abs(strain[0]+emax) < 1.e-7):
                strain.pop(0)
                energy.pop(0)
            if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                strain.pop()
                energy.pop()

    fD.close()
    fE.close()


    #--- Plotting -----------------------------------------------------------------------------
    if (os.path.exists('Grace.par') == False):
        os.system("cp -f $EXCITINGSCRIPTS/Grace.par .")

    Gf    = open('Grace.par', 'r')
    Glines= Gf.readlines()
    Gf.close()

    TMP = []
    for k in range(221, 323):
        TMP.append(Glines[k])

    TMP.insert(46,'    subtitle "Plot for '+ Dstn +' deformation, n = Order of polynomial fit"\n')

    GdE = open(Dstn+'_d1E.par', 'w')
    for k in range(len(TMP)):
        print >>GdE, TMP[k],
    GdE.close()

    os.system('xmgrace '+ Dstn + '_d1E.dat -param '  + \
                          Dstn + '_d1E.par -saveall '+ \
                          Dstn + '_d1E.agr &')

    TMP = []
    for k in range(154, 162):
        TMP.append(Glines[k])

    for k in range(164, 219):
        TMP.append(Glines[k])

    TMP.insert(63,'    s2 legend  " n = 1 "\n')
    TMP.insert(55,'    s1 legend  " n = 3 "\n')
    TMP.insert(47,'    s0 legend  " n = 5 "\n')
    TMP.insert(10,'    subtitle "Plot for '+ Dstn +' deformation, n = Order of polynomial fit"\n')

    CVe = open(Dstn+'_CVe.par', 'w')
    for k in range(len(TMP)):
        print >>CVe, TMP[k],
    CVe.close()

    #os.system('xmgrace '+ Dstn +'_CVe.dat -param '+ Dstn +'_CVe.par -saveall '+ Dstn +'_CVe.agr &')

os.chdir('../')

#--- Writing the "STRESS.IN" file -----------------------------------------------------------------

fri = open('STRESS.IN', 'w')
for i in range(1, SCs+1):
    if (i<10):
        Dstn = 'Dst0'+str(i)
    else:
        Dstn = 'Dst' +str(i)
    print >>fri, Dstn+'    eta_max    Fit_order'

fri.close()
#--------------------------------------------------------------------------------------------------
os.system('rm -f Energy-vs-Strain/Grace.par')
