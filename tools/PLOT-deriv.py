#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   sys   import stdin
from   math  import sqrt
from   math  import factorial
from   pylab import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ptk 
import pylab             as pyl
import numpy
import sys
import os  
import math 
import numpy
from   scipy import *
from   scipy.optimize import leastsq

#-------------------------------------------------------------------------------

def sortstrain(s,e):
    ss=[] ; ee=[] ; ww=[]
    for i in range(len(s)): ww.append(s[i])
    ww.sort()
    for i in range(len(s)):
        ss.append(s[s.index(ww[i])])
        ee.append(e[s.index(ww[i])])
    return ss, ee
    
#-------------------------------------------------------------------------------

def shell_value(variable,vlist,default):
    v = default
    e = False
    for i in range(len(vlist)):
        if ( vlist[i] == variable ): v = os.environ[variable] ; e = True ; break
    return v, e
    
#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

if (narg < 1): 
    print "\nIncorrect number of arguments. **Usage**:\n\n",
    print "PLOT-deriv.py FILENAME\n"
    sys.exit()

inputfile = str(sys.argv[1])
lexsev    = os.path.exists(inputfile)
if (not(lexsev)):   sys.exit("ERROR: file "+inputfile+" not found!\n")

input_energy = open(inputfile,"r")

energy = [] ; volume = []

while True:
    line = input_energy.readline().strip()
    if len(line) == 0: break
    energy.append(float(line.split()[1]))
    volume.append(float(line.split()[0]))

nv = len(volume)
if (nv < 2): sys.exit("\nERROR: Too few x values ("+str(nv)+")!\n")

#-------------------------------------------------------------------------------

volume, energy = sortstrain(volume,energy)

#-------------------------------------------------------------------------------

input_energy.close()

#*******************************************************************************

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])

#-------------------------------------------------------------------------------

xlabel = u'x'
ylabel = r'y'

#-------------------------------------------------------------------------------

pvol = [] ; ppre = []
for i in range(len(volume)-1): 
    pppp=(energy[i+1]-energy[i])/(volume[i+1]-volume[i])
    ppre.append(pppp)
    pvol.append(volume[i]+(volume[i+1]-volume[i])/2)

#-------------------------------------------------------------------------------

fontlabel=20
fonttick=16

params = {'ytick.minor.size': 6,
          'xtick.major.pad': 8,
          'ytick.major.pad': 4,
          'patch.linewidth': 2.,
          'axes.linewidth': 2.,
          'lines.linewidth': 1.8,
          'lines.markersize': 8.0,
          'axes.formatter.limits': (-4, 6)}

plt.rcParams.update(params)
plt.subplots_adjust(left=0.21, right=0.93,
                    bottom=0.18, top=0.88,
                    wspace=None, hspace=None)
                           
yfmt = ptk.ScalarFormatter(useOffset=True,useMathText=True)                           
                           
figure = plt.figure(1, figsize=(8,5.5))  
ax     = figure.add_subplot(111)
ax.text(0.5,-0.13,xlabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center')
ax.text(-0.19,0.5,ylabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=90)
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(6)
    line.set_markeredgewidth(2)
plt.xticks(size=fonttick)
plt.yticks(size=fonttick)
pyl.grid(True)

plt.plot(pvol,ppre,'r-')
plt.plot(pvol,ppre,'go',label='finite differences')
#plt.legend(loc=4,borderaxespad=.8,numpoints=1)
plt.legend(bbox_to_anchor=(1.0, 1.03), loc=4, borderaxespad=0.)

ymax  = max(ppre)
ymin  = min(ppre)
dxx   = abs(max(volume)-min(volume))/18
dyy   = abs(ymax-ymin)/18
ax.yaxis.set_major_formatter(yfmt)
ax.set_xlim(min(volume)-dxx,max(volume)+dxx)
ax.set_ylim(ymin-dyy,ymax+dyy)

ax.xaxis.set_major_locator(MaxNLocator(7))

plt.savefig('PLOT.ps', orientation='portrait',format='eps')
plt.savefig('PLOT.png',orientation='portrait',format='png',dpi=dpipng)

#-------------------------------------------------------------------------------





