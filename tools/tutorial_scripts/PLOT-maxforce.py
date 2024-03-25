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

xlabel  = u'Optimization steps'
ylabel  = r'Maximum force [Ha/Bohr]'

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
    print "PLOT-maxforce.py DIRECTORYNAME\n"
    sys.exit()

label = str(sys.argv[1])

#-------------------------------------------------------------------------------

inpf = label+'/INFO.OUT'
     
if (str(os.path.exists(inpf))=='False'): 
    sys.exit("\nERROR: file "+inpf+" not found!\n")
    
#-------------------------------------------------------------------------------

os.system("grep \"Maximum force m\" "+str(inpf)+" > tempfile")
input_file = open("tempfile","r")

if ( os.path.getsize('./tempfile') == 0 ):
    print "\nEither data file not (yet) ready for visualization",
    print "\nor maximum force target reached already at the initial configuration.\n"
    sys.exit() 

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
a = [] ; b = []

iter=0
lgoal=0

while True:
    line = input_file.readline().strip().replace(")", "")
    if len(line) == 0: break
    iter+=1
    y.append(float(line.split()[5]))
    x.append(float(iter-1))
    if (iter == 1):
        lgoal=1
        goal=float(line.split()[7])
        b.append(goal)
        a.append(float(iter-1))

if (abs(y[-1]) < 1.e-8 ): y[-1] = 1.e-8
        
if (lgoal < 1): 
    os.system("rm -f tempfile")
    print "\nData not (yet) available for visualization.\n\n",
    sys.exit() 
   
b.append(goal)
a.append(float(iter-1))

os.system("rm -f tempfile")

xmin = 0-(iter-1)/20. ; xmax = (iter-1)+(iter-1)/20.

if (iter == 1): 
    xmin = -1 ; xmax = 1

#-------------------------------------------------------------------------------

plt.plot(a,b,'b-',label='goal')
plt.plot(x,y,'ro-',label='run')

#-------------------------------------------------------------------------------

plt.legend(borderaxespad=.8)

ax.yaxis.set_major_formatter(yfmt)
ax.set_yscale('log')
ax.set_xlim(xmin,xmax)

ax.set_axisbelow(True) 

plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------





