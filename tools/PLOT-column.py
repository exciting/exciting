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

def sortstrain(s,e):
    ss=[]
    ee=[]
    ww=[]
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

def leggi(filin,c1,c2):
    onlyone = False
    f = open(filin,"r")
    x = [] ; y = []
    cmax=max(c1,c2)
    ix = 0
    lines = f.readlines()
    for ilines in range(len(lines)):
        line = lines[ilines].strip() 
        if (len(line) != 0):
            if (len(line) < cmax): 
                sys.exit("\n ERROR: number of columns is smaller than "+str(cmax)+"!\n")
            if (c1==c2): 
                onlyone = True
                x.append(int(ix))
                y.append(float(line.split()[c2-1])) 
            else:
                x.append(float(line.split()[c1-1]))
                y.append(float(line.split()[c2-1])) 
            ix = ix+1
    f.close()
    return x,y,onlyone

#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])

#-------------------------------------------------------------------------------

fontlabel=20
fonttick=16
fonttext=14

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

fig  = matplotlib.pyplot.figure(1, figsize=(8,5.5)) 

ax   = fig.add_subplot(111)

#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

if (narg<1): 
    print "\nIncorrect number of arguments. **Usage**:\n\n",
    print "PLOT-column.py datafile [COLUMN1 COLUMN2 POINTSTYLE YMIN  [YMAX]]\n"
    sys.exit()

filin = str(sys.argv[1])

column1 = 1
column2 = 2

if (narg==2): 
    column1 = int(sys.argv[2])
    column2 = int(sys.argv[2])

if (narg>2):
    column1 = int(sys.argv[2])
    column2 = int(sys.argv[3])

pointstyle='o-'
if (narg>3): pointstyle=str(sys.argv[4])

ylimits = []
for i in range(5,len(sys.argv)): ylimits.append(float(sys.argv[i]))

#-------------------------------------------------------------------------------

xmin = 1.e30 ; xmax = -1.e30
ymin = 1.e30 ; ymax = -1.e30

x = [] ; y = []
x, y, onlyone = leggi(filin,column1,column2)

plt.plot(x,y,'r'+pointstyle,label="file")
xmin=min(min(x),xmin) ; xmax=max(max(x),xmax)
ymin=min(min(y),ymin) ; ymax=max(max(y),ymax)
    
#-------------------------------------------------------------------------------
# set defauls parameters for the plot
 
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(6)
    line.set_markeredgewidth(2)      

plt.xticks(size=fonttick)
plt.yticks(size=fonttick)
pyl.grid(True)

#-------------------------------------------------------------------------------

plt.legend(loc=1,borderaxespad=.8,numpoints=1)

ax.yaxis.set_major_formatter(yfmt)

dxx = (xmax-xmin)/18 ; dyy = (ymax-ymin)/15
xmin = xmin-dxx ; xmax = xmax+dxx
ymin = ymin-dyy ; ymax = ymax+dyy

ylimits = []
for i in range(5,len(sys.argv)): ylimits.append(float(sys.argv[i]))

if (len(ylimits) == 1): ymin = float(ylimits[0])
if (len(ylimits) > 1): ymin = float(ylimits[0]); ymax = float(ylimits[1]) 

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

ax.xaxis.set_major_locator(MaxNLocator(7))

ax.set_axisbelow(True) 

ax.text(1,1.05,'file: '+filin,size=fontlabel,
        transform=ax.transAxes,ha='right',va='center',rotation=0)
        
plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------





