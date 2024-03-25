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

def shell_value(variable,vlist,default):
    v = default
    e = False
    for i in range(len(vlist)):
        if ( vlist[i] == variable ): v = os.environ[variable] ; e = True ; break
    return v, e

#-------------------------------------------------------------------------------

def readtext(filin,text,itext):
    os.system("grep \""+text+"\" "+filin+" | tail -n1 > tempfile")
    f = open("tempfile","r")
    e = f.readline().strip().split()[itext]
    os.system("rm -f tempfile")
    return e
        
#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])
   
#-------------------------------------------------------------------------------

#print "\n**Usage**:    PLOT-genenergy.py [DIRECTORYROOT]\n"

#-------------------------------------------------------------------------------

directoryroot = 'rundir-'
if (len(sys.argv) > 1): directoryroot = str(sys.argv[1])
list_dir = sorted(glob.glob(directoryroot+"*"))

#-------------------------------------------------------------------------------

#if ( not(os.path.exists(current+'/'+list_dir[0]+'/input.xml')) ): 
#    sys.exit("\n ERROR: Input file "+current+"/"+list_dir[0]+"/input.xml not found!\n")

#-------------------------------------------------------------------------------

xx = [] ; yy = []

text  = "Total energy                               :"
itext = 3

for idir in range(len(list_dir)):
    fileinp = list_dir[idir]+'/INFO.OUT'
    if ( not(os.path.exists(fileinp)) ): break
    yy.append(float(readtext(fileinp,text,itext)))
    xx.append(float(idir+1))
          
ylabel = r'Energy'
xlabel = r'run index'

#-------------------------------------------------------------------------------
# manipulate data for a better plot

xmin = min(xx) ; xmax = max(xx)
ymin = min(yy) ; ymax = max(yy)

dxx  = abs(xmax-xmin)/18
dyy  = abs(ymax-ymin)/15

xmin = xmin-dxx ; xmax = xmax+dxx
ymin = ymin-dyy ; ymax = ymax+dyy

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

plt.plot(xx,yy,'ro--')

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

#if (len(sys.argv) > 4): ymin = float(sys.argv[4])
#if (len(sys.argv) > 5): ymax = float(sys.argv[5])

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

ax.xaxis.set_major_locator(MaxNLocator(7))

ax.set_axisbelow(True) 

#plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------





