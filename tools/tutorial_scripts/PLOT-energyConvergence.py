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

xlabel  = u'Iteration number $k$'
ylabel  = r'|E(k)-E(k$_{last}$)| [Ha]'

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
    print "PLOT-energyConvergence.py DIRECTORYNAME\n"
    sys.exit()

label = str(sys.argv[1])

#-------------------------------------------------------------------------------

rsmd_file = label+'/TOTENERGY.OUT'
info_file = label+'/INFO.OUT'
lrsmd = os.path.exists(rsmd_file)
linfo = os.path.exists(info_file)
lonlyinfo=False

if (lrsmd):
    icol = 0
    input_file = open(rsmd_file,"r")
else:
    sys.exit("\nWARNING: file "+info_file+" not (yet) found!\n")
    
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

plt.rcParams.update({'mathtext.default':'regular'})

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

ene = []

while True:
    line = input_file.readline().strip().replace(")", "") 
    if len(line) != 0:
       ene.append(float(line.split()[0]))
    else: 
       break

x = []
y = []

for i in range(len(ene)-1):
    x.append(float(i+1))
    y.append(abs(float(ene[i]-ene[-1])))
       
plt.plot(x,y,'b-')
plt.plot(x,y,'go',label='calculated')
plt.legend(borderaxespad=.8,numpoints=1,framealpha=0.9,fancybox=True)

dx = x[-1]- x[0]
xmin = 1-dx/20. ; xmax = x[-1]+dx/20.

#-------------------------------------------------------------------------------
   
ax.yaxis.set_major_formatter(yfmt)
ax.set_yscale('log')
ax.set_xlim(xmin,xmax)

ax.set_axisbelow(True) 

plt.savefig('PLOT.png',format='png',dpi=300, bbox_inches='tight')
plt.savefig('PLOT.eps',format='eps',bbox_inches='tight')

#plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
#plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------





