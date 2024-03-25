#!/usr/bin/env python
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%%% -------------------------------- ElaStic@exciting-result -------------------------------- %%%#
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
# python ElaStic@exciting-result.py
# 
# EXPLANATION:
# 
#__________________________________________________________________________________________________

import sys
import os

#%!%--- Checking the INFO_ElaStic exist ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
if (os.path.exists('INFO_ElaStic') == False):
    sys.exit('\n     ... Oops ERROR: Where is the "INFO_ElaStic" file !?!?!?    \n')
#--------------------------------------------------------------------------------------------------

#%!%--- Reading the INFO_ElaStic file ---#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
INFO=open('INFO_ElaStic', 'r')

l1  = INFO.readline()
ordr= int(l1.split()[-1])

if (ordr != 2 and ordr != 3):
    sys.exit('\n     ... Oops ERROR: The order of the elastic constant is NOT clear !?!?!?'\
             '\n                     Something is WRONG in "INFO_ElaStic" file.\n')

l2  = INFO.readline()
mthd= l2.split()[-1]

if (mthd != 'Energy'):
    sys.exit('\n     ... Oops ERROR: The method of the calculation is NOT clear !?!?!?'\
             '\n                     Something is WRONG in "INFO_ElaStic" file.\n')

INFO.close()
#--------------------------------------------------------------------------------------------------

#%!%--- Checking the ElaStic_???.in exist ---#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (ordr == 2):
    if (os.path.exists('ElaStic_2nd.in') == False):
        sys.exit('\n     ... Oops ERROR: Where is the "ElaStic_2nd.in" file !?!?!?\n')

if (ordr == 3):
    if (os.path.exists('ElaStic_3rd.in') == False):
        sys.exit('\n     ... Oops ERROR: Where is the "ElaStic_3rd.in" file !?!?!?\n')

#--------------------------------------------------------------------------------------------------

if (mthd == 'Energy' and ordr == 2): os.system('ElaStic_Result_Energy_2nd.py')
if (mthd == 'Energy' and ordr == 3): os.system('ElaStic_Result_Energy_3rd.py')
