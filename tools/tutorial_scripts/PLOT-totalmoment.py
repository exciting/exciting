#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   lxml  import etree
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

def leggi(filin,a):
    os.system("grep \""+a+"\" "+str(filin)+" | tail -n1 > tempfile") 
    ifile = open("tempfile","r")
    y = ifile.readline().strip().split()[3]
    ifile.close()
    os.system("rm -f tempfile")
    return y
    
#-------------------------------------------------------------------------------

def readstrain(dr):
    ifx = open("energy-vs-volume","r")   
    x = []
    while True:
        line = ifx.readline().strip()
        if len(line) == 0: break
        x.append(float(line.split()[0]))
    ifx.close()
    return x
    
#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])
   
#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

print "\n**Usage**:    PLOT-totalmoment.py [DIRECTORYROOT]\n"

#-------------------------------------------------------------------------------

directoryroot = 'rundir-'
if (len(sys.argv) > 1): directoryroot = str(sys.argv[1])
list_dir = sorted(glob.glob(directoryroot+"*"))

#-------------------------------------------------------------------------------

x = [] ; y = []

for idir in range(len(list_dir)):
    fileinp = list_dir[idir]+'/INFO.OUT'
    r = leggi(fileinp,"total moment    ")
    y.append(float(r))
    
x = readstrain(directoryroot)
        
ylabel  = r'Total moment [$\mu_B$]'
xlabel  = r'Volume'

#-------------------------------------------------------------------------------
# manipulate data for a better plot

xmin = min(x)
xmax = max(x)

ymin = min(y)
ymax = max(y)

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

plt.plot(x,y,'ro-')

#-------------------------------------------------------------------------------

#plt.legend(loc=9,borderaxespad=.8,numpoints=1)
#plt.legend(bbox_to_anchor=(1.03, 1), loc=2, borderaxespad=0., numpoints=1)

ax.yaxis.set_major_formatter(yfmt)

#ylimits = []
#for i in range(1,len(sys.argv)): ylimits.append(float(sys.argv[i]))

#if (len(ylimits) == 1): ymin = float(ylimits[0])
#if (len(ylimits) > 1): ymin = float(ylimits[0]); ymax = float(ylimits[1]) 

if (abs(ymax-ymin) < 0.000000001): 
    ymax=ymax+0.1
    ymin=ymin-0.1

if (len(sys.argv) > 4): ymin = float(sys.argv[4])
if (len(sys.argv) > 5): ymax = float(sys.argv[5])

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

ax.xaxis.set_major_locator(MaxNLocator(7))

ax.set_axisbelow(True) 

plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

output_totmoment = open('total-moment',"w")

pmt='%16.10f' 

for i in range(len(x)): 
    print >> output_totmoment, pmt%(x[i]), pmt%(y[i])

output_totmoment.close()

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------





