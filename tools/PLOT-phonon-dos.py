#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   lxml  import etree as et
from   pylab import *
import sys
import matplotlib         as mpl
import matplotlib.pyplot  as plt
import matplotlib.ticker  as ptk 
import pylab              as pyl
import numpy
import os

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

#-------------------------------------------------------------------------------
# Read arguments

narg   = len(sys.argv)-1
plabel = "|"
mcolor = "w"

if (narg > 0): xxmin  = float(sys.argv[1])
if (narg > 1): xxmax  = float(sys.argv[2])
if (narg > 2): 
    plabel = str(sys.argv[3])
    mcolor = "b" 
   
#-------------------------------------------------------------------------------
# Read data

ene=[]
dos=[]
for line in open("PHDOS.OUT"):
    i_line=line.split()
    if len(i_line):
       eneval=float(i_line[0])*2.194746313705e5
       dosval=float(i_line[1])/2.194746313705e5
       if ( eneval > 0.):
          ene.append(eneval)
          dos.append(dosval)

#-------------------------------------------------------------------------------
# Settings for the plot

fontlabel=22
fontlegend=16
fonttick=16
fonttickx=16

params = {'font.family': 'sans-serif',
          'ytick.minor.size': 6,
          'xtick.major.pad': 8,
          'ytick.major.pad': 4,
          'patch.linewidth': 2.,
          'axes.linewidth': 2.5,
          'lines.linewidth': 2.0,
          'lines.markersize': 5.0,
          'axes.formatter.limits': (-5, 6)}

plt.rcParams.update(params)

plt.subplots_adjust(left=0.20, right=0.93,
                    bottom=0.28, top=0.92,
                    wspace=None, hspace=None)
                    
yfmt = ptk.ScalarFormatter(useOffset=True,useMathText=True)
fig  = matplotlib.pyplot.figure(1, figsize=(8,5.5)) 
ax1  = fig.add_subplot(111)

ax1.xaxis.set_label_position('bottom')

ax1.text(-0.18,0.5,'Phonon DOS [states/cm$\mathregular{^{-1}}$]',size=fontlabel,
        transform=ax1.transAxes,ha='center',va='center',rotation=90)
ax1.text(0.5,-0.16,'Frequency [cm$\mathregular{^{-1}}$]',size=fontlabel,
        transform=ax1.transAxes,ha='center',va='center',rotation=0)

#-------------------------------------------------------------------------------
# Tick size

for line in ax1.get_xticklines() + ax1.get_yticklines(): line.set_markersize(9)
for line in ax1.get_yticklines(): line.set_markeredgewidth(2)
for line in ax1.get_xticklines(): line.set_markeredgewidth(2)

#-------------------------------------------------------------------------------
# Plot

ax1.plot(ene,dos,'r-',label="",color="firebrick")
plt.fill_between(ene,dos, color="darksalmon")
		 
if ( narg > 0): plt.xlim(xmin=xxmin)
if ( narg > 1): plt.xlim(xmax=xxmax)

plt.xticks(size=fonttickx)
plt.yticks(size=fonttick)

ax1.text(1,1.05,plabel,size=fontlabel, color=mcolor,
             transform=ax1.transAxes,ha='right',va='center',rotation=0)

fig.subplots_adjust(left=None, bottom=None, right=None, wspace=None, hspace=None)

fig.savefig('PLOT.png',format='png',bbox_inches='tight',dpi=300)
fig.savefig('PLOT.pdf',format='pdf',bbox_inches=0)

sys.exit()    

#-------------------------------------------------------------------------------