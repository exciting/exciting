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

if (os.path.exists('exciting')):         ylabel = r'Energy [Ha]'
if (os.path.exists('quantum-espresso')): ylabel = r'Energy [Ry]'
if (os.path.exists('vasp')):             ylabel = r'Energy [Ry]'

inpf = ' '

if (os.path.exists('energy-vs-strain')): 
     inpf   = 'energy-vs-strain'
     xlabel = r'Lagrangian strain'

if (os.path.exists('energy-vs-displacement')): 
    inpf   = 'energy-vs-displacement'
    xlabel = r'Displacement $u$ [alat]'

if (os.path.exists('energy-vs-volume')): 
    inpf   = 'energy-vs-volume'
    xlabel = u'Volume [Bohr\u00B3]'
     
if (os.path.exists('energy-vs-ngrid')): 
    inpf   = 'energy-vs-ngridk'
    xlabel = r'ngridk'
     
if (os.path.exists('energy-vs-alat')): 
    inpf   = 'energy-vs-alat'
    xlabel = r'Lattice parameter [Bohr]'
     
if (str(os.path.exists(inpf))=='False'): 
    sys.exit("\nERROR: file "+inpf+" not found!\n")
    
if (os.path.exists('xlabel')): 
    label_file = open('xlabel',"r")
    xlabel = str(label_file.readline())
    label_file.close()

#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])

#-------------------------------------------------------------------------------

input_file = open(inpf,"r")

x = [] ; y = []

while True:
    line = input_file.readline()
    line = line.strip()
    if len(line) == 0: break
    y.append(float(line.split()[1]))
    x.append(float(line.split()[0]))

xx,yy=sortstrain(x,y)

#-------------------------------------------------------------------------------

if (os.path.exists('rundir-oo')):
    infty_file = open("./rundir-oo/TOTENERGY.OUT","r")
    energy_zero = float(infty_file.readlines()[-2].strip())
    infty_file.close()
    yy.pop(-1)
    xx.pop(-1)
    norm_file = open("normalized-energy","w")
    for i in range(len(yy)): 
        yy[i] = round(yy[i]-energy_zero,8)
        print >>norm_file, xx[i], yy[i]      
    norm_file.close()
        
#-------------------------------------------------------------------------------
# manipulate data for a better plot

rmin  = min(yy)
srmin = ""

if (not(os.path.exists('rundir-oo'))):
    srmin = u'\u2013 '+str(abs(rmin))
    if (rmin > 0): srmin = u'+ '+str(rmin)
    for i in range(len(yy)): yy[i]=(yy[i]-rmin)

delta = (max(yy)-min(yy))
dxx   = abs(max(xx)-min(xx))/18
dyy   = abs(max(yy)-min(yy))/17

xmin = min(xx)-dxx ; xmax = max(xx)+dxx
ymin = min(yy)-dyy ; ymax = max(yy)+dyy

if (len(xx) == 1): 
    xmin = min(xx)-dxx-1 ; xmax = max(xx)+dxx+1
    ymin = min(yy)-dyy ; ymax = max(yy)+dyy

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
plt.subplots_adjust(left=0.21, right=0.93,
                    bottom=0.18, top=0.88,
                    wspace=None, hspace=None)
                    
yfmt = ptk.ScalarFormatter(useOffset=True,useMathText=True)

fig  = matplotlib.pyplot.figure(1, figsize=(8,5.5)) 

ax   = fig.add_subplot(111)

if (os.path.exists('rundir-oo')):       
    ax.text(0.5,-0.15,xlabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=0)
else:
    ax.text(0.5,-0.14,xlabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=0)
    ax.text(0.11,1.03,srmin,size=fonttext,
        transform=ax.transAxes,ha='left',va='center',rotation=0)

ax.text(-0.165,0.5,ylabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=90)
       
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(6)
    line.set_markeredgewidth(2)      

plt.xticks(size=fonttick)
plt.yticks(size=fonttick)
pyl.grid(True)

plt.plot(xx,yy,'r-')
plt.plot(xx,yy,'go',label='calculated')

#-------------------------------------------------------------------------------

plt.legend(loc=9,borderaxespad=.8,numpoints=1)

ax.yaxis.set_major_formatter(yfmt)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

ax.xaxis.set_major_locator(MaxNLocator(7))

ax.set_axisbelow(True) 

plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------





