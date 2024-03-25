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

def leggi(filin):
    onlyone = False
    f = open(filin,"r")
    x = [] ; y = []
    ilines = 0
    ix = 0
    lines = f.readlines()
    while ( ilines < len(lines) ):
        line = lines[ilines].strip() 
        ilines = ilines+1
        if (len(line) != 0 and line[0][:1] != "#"): 
            if (len(line.split()) == 1):
	        onlyone = True
                x.append(int(ix))
                y.append(float(line.split()[0])) 
            else:
                x.append(float(line.split()[0]))
                y.append(float(line.split()[1])) 
            ix = ix+1
            
    #ymax=max(y)
    #print ymax
    #z = []
    #for i in range(len(y)):
    #    z.append(float(y[i]/ymax))
            
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
plt.subplots_adjust(left=0.21, right=0.78,
                    bottom=0.18, top=0.88,
                    wspace=None, hspace=None)
                    
yfmt = ptk.ScalarFormatter(useOffset=True,useMathText=True)

fig  = matplotlib.pyplot.figure(1, figsize=(8,5.5)) 

ax   = fig.add_subplot(111)

#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

if (narg<1 or narg>2): 
    print "\nIncorrect number of arguments. **Usage**:\n\n",
    print "PLOT-diff.py FILE1 FILE2\n"
    sys.exit()

#-------------------------------------------------------------------------------
# set defauls parameters for the plot
 
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(6)
    line.set_markeredgewidth(2)      

plt.xticks(size=fonttick)
plt.yticks(size=fonttick)
pyl.grid(True)

#-------------------------------------------------------------------------------

xmin = 1.e30 ; xmax = -1.e30
zmin = 1.e30 ; zmax = -1.e30
 
superx = [] ; supery = [] 
filin1 = str(sys.argv[1])
filin2 = str(sys.argv[2])

x1 = [] ; y1 = []
x1, y1, onlyone = leggi(filin1)  
    
x2 = [] ; y2 = []
x2, y2, onlyone = leggi(filin2)  

z = [] 
for i in range(len(y1)): 
    z.append(float(y2[i]-y1[i]))    
    
asym = "-o"
if (len(x1) > 50): asym = "-" 
    #for i in range(1,len(y)): y[i] = y[i]-y[0]
    #y[0] = 0.0
plt.plot(x1,z,asym)#"file"+str(ifile+1))

xmin=min(min(x1),xmin) ; xmax=max(max(x1),xmax)
zmin=min(min(z), zmin) ; zmax=max(max(z), zmax)
   
#-------------------------------------------------------------------------------

#xmax=0.2

#plt.legend(loc=1,borderaxespad=.8,numpoints=1)
#plt.legend(bbox_to_anchor=(1.03, 1), loc=2, borderaxespad=0., numpoints=1)

ax.yaxis.set_major_formatter(yfmt)

dxx = (xmax-xmin)/18 ; dzz = (zmax-zmin)/15
xmin = xmin-dxx ; xmax = xmax+dxx
zmin = zmin-dzz ; zmax = zmax+dzz

ax.set_xlim(xmin,xmax)
ax.set_ylim(zmin,zmax)

ax.xaxis.set_major_locator(MaxNLocator(7))

ax.set_axisbelow(True) 

#ax.text(1,1.05,'file: '+filin,size=fontlabel,
#        transform=ax.transAxes,ha='right',va='center',rotation=0)
        
#plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------





