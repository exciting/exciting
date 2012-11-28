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
import glob
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

def leggi(filin):
    f = open(filin,"r")
    x = [] ; y = []
    lines = f.readlines()
    y = lines[-1].strip().split()[0:3]
    x = lines[-2].strip().split()[0:3]
    f.close()
    return x,y

#-------------------------------------------------------------------------------

def leggishort(filin):
    f = open(filin,"r")
    x = float(f.readline().strip().split()[0])
    f.close()
    return x

#-------------------------------------------------------------------------------

list_dir = glob.glob('rundir-*')

xx = [] ; x0 = []
y1 = [] ; y2 = [] ; y3 = [] 
z1 = [] ; z2 = [] ; z3 = [] 

for idir in range(len(list_dir)):
    r1 = [] ; r2 = []    
    c1 = [] ; c2 = []
    fileref = list_dir[idir]+'/GEOMETRY.OUT'
    fileopt = list_dir[idir]+'/GEOMETRY_OPT.OUT'
    if (str(os.path.exists(fileopt))=='False'): fileopt = fileref 
    r1,r2 = leggi(fileref)
    c1,c2 = leggi(fileopt)
    y1.append(float(c2[0])-float(c1[0]))
    y2.append(float(c2[1])-float(c1[1]))
    y3.append(float(c2[2])-float(c1[2]))
    xx.append(leggishort('strain-'+list_dir[idir][-2:]))
    if (idir == 0 or idir == len(list_dir)-1):
        z1.append(float(r2[0])-float(r1[0]))
        z2.append(float(r2[1])-float(r1[1]))
        z3.append(float(r2[2])-float(r1[2]))
        x0.append(leggishort('strain-'+list_dir[idir][-2:]))
    
ylabel  = r'Relative coordinate [crystal]'
xlabel  = r'Lagrangian strain'

#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])

#-------------------------------------------------------------------------------
# manipulate data for a better plot

xmin = min(xx+x0)
xmax = max(xx+x0)

ymin = min(y1+y2+y3+z1+z2+z3)
ymax = max(y1+y2+y3+z1+z2+z3)

dxx  = abs(xmax-xmin)/18
dyy  = abs(ymax-ymin)/15

xmin = xmin-dxx
xmax = xmax+dxx

ymin = ymin-dyy
ymax = ymax+dyy

#-------------------------------------------------------------------------------
# set defauls parameters for the plot

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
plt.subplots_adjust(left=0.22, right=0.78,
                    bottom=0.18, top=0.88,
                    wspace=None, hspace=None)
                    
yfmt = ptk.ScalarFormatter(useOffset=True,useMathText=True)

fig  = matplotlib.pyplot.figure(1, figsize=(8,5.5)) 

ax   = fig.add_subplot(111)

ax.text(0.5,-0.13,xlabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=0)
ax.text(-0.23,0.5,ylabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=90)
 
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(6)
    line.set_markeredgewidth(2)      

plt.xticks(size=fonttick)
plt.yticks(size=fonttick)
pyl.grid(True)

plt.plot(xx,y1,'ro--',label=u'$\Delta$1')
plt.plot(xx,y2,'bo--',label=u'$\Delta$2')
plt.plot(xx,y3,'go--',label=u'$\Delta$3')

plt.plot(x0,z1,'r-',label=u'$\Delta$1$_{ref}$')
plt.plot(x0,z2,'b-',label=u'$\Delta$2$_{ref}$')
plt.plot(x0,z3,'g-',label=u'$\Delta$3$_{ref}$')

plt.plot(x0,z3,'g-')
plt.plot(x0,z2,'b-')
plt.plot(x0,z1,'r-')

plt.plot(xx,y3,'go--')
plt.plot(xx,y2,'bo--')
plt.plot(xx,y1,'ro--')

#-------------------------------------------------------------------------------

#plt.legend(loc=9,borderaxespad=.8,numpoints=1)
plt.legend(bbox_to_anchor=(1.03, 1), loc=2, borderaxespad=0.)

ax.yaxis.set_major_formatter(yfmt)

ylimits = []
for i in range(1,len(sys.argv)): ylimits.append(float(sys.argv[i]))

if (len(ylimits) == 1): ymin = float(ylimits[0])
if (len(ylimits) > 1): ymin = float(ylimits[0]); ymax = float(ylimits[1]) 

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

ax.xaxis.set_major_locator(MaxNLocator(7))

ax.set_axisbelow(True) 

plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------





