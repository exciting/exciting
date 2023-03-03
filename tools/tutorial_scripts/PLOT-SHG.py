#!/usr/bin/python2
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   sys   import stdin
from   math  import sqrt
from   math  import factorial
from   pylab import *
import matplotlib.transforms as ptf
import matplotlib.ticker     as ptk 
import matplotlib.pyplot     as plt
import pylab                 as pyl
import numpy
import sys
import os

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

#-------------------------------------------------------------------------------

def shell_value(variable,vlist,default):
    v = default
    e = False
    for i in range(len(vlist)):
        if ( vlist[i] == variable ): v = os.environ[variable] ; e = True ; break
    return v, e

#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])

#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

if (narg<1): 
    print "\n Incorrect number of arguments. **Usage**:\n\n",
    print " PLOT-SHG.py datafile\n"
    sys.exit()

filin = str(sys.argv[1])

options = " -cy 2 3 4 -ll Real Imaginary Modulus -lx 'Energy [eV]' -g -rc"
options = options + " -ly '$\chi^{(2)}(-2\omega,\omega,\omega)$ [$10^{-7}$ esu]'"
os.system("PLOT-files.py -f "+filin+" "+filin+" "+filin+options)

#-------------------------------------------------------------------------------





