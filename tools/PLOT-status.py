#!/usr/bin/python
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

#-------------------------------------------------------------------------------

def shell_value(variable,vlist,default):
    v = default
    e = False
    for i in range(len(vlist)):
        if ( vlist[i] == variable ): v = os.environ[variable] ; e = True ; break
    return v, e

#-------------------------------------------------------------------------------

xlabel  = u'Iteration number'
ylabel  = r'RMSD $V_{SCF}$'

#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])

#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

if (narg < 1): 
    print "\nIncorrect number of arguments. **Usage**:\n\n",
    print "PLOT-status.py DIRECTORYNAME\n"
    sys.exit()

label = str(sys.argv[1])

#-------------------------------------------------------------------------------

inpf = current+"/"+rlabel+label+'/RMSDVEFF.OUT'
icol = 0

if (str(os.path.exists(inpf))=='False'): 
    inpf    = current+"/"+rlabel+label+'/INFO.OUT'
    if (str(os.path.exists(inpf))=='False'): 
        sys.exit("\nERROR: file "+inpf+" not found!\n")
    icol = 7
    os.system("grep \"RMS change\" "+str(inpf)+" > tempfile")
    input_file = open("tempfile","r")
else:
    input_file = open(inpf,"r")

#-------------------------------------------------------------------------------
# set defauls parameters for the plot

fontlabel=20
fonttick=16

params = {'ytick.minor.size': 6,
          'xtick.major.pad': 8,
          'ytick.major.pad': 4,
          'patch.linewidth': 2.,
          'axes.linewidth': 2.,
          'lines.linewidth': 1.8,
          'lines.markersize': 8.0,
          'axes.formatter.limits': (-5, 6)}

plt.rcParams.update(params)

plt.subplots_adjust(left=0.20, right=0.93,
                    bottom=0.18, top=0.88,
                    wspace=None, hspace=None)
                    
yfmt = ptk.ScalarFormatter(useOffset=True,useMathText=True)

fig  = matplotlib.pyplot.figure(1, figsize=(8,5.5)) 

ax   = fig.add_subplot(111)

ax.text(0.5,-0.13,xlabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=0)
ax.text(-0.16,0.5,ylabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=90)

for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(6)
    line.set_markeredgewidth(2)
       
plt.xticks(size=fonttick)
plt.yticks(size=fonttick)
pyl.grid(True)

#-------------------------------------------------------------------------------

x = [] ; y = []

iter=1
alast=0
alege=0

while True:
    line = input_file.readline().strip().replace(")", "") 
    if len(line) == 0:
       if alast == 1: break
       plt.plot(x,y,'b-')
       plt.plot(x,y,'go',label='calculated')
       if alege == 0: plt.legend(borderaxespad=.8,numpoints=1)
       x = [] ; y = []
       alast = 1
       alege = 1
    else:
       alast = 0
       iter  = iter+1
       y.append(float(line.split()[icol]))
       x.append(float(iter))
 
if (str(os.path.exists(inpf))=='False'): os.system("rm tempfile")

xmin = 2-iter/20. ; xmax = iter+iter/20.

#-------------------------------------------------------------------------------
   
ax.yaxis.set_major_formatter(yfmt)
ax.set_yscale('log')
ax.set_xlim(xmin,xmax)

ax.set_axisbelow(True) 

plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------





