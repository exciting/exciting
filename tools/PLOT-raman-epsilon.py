#!/usr/bin/env python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   lxml  import etree as et
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as pyl
import os
import numpy as np
from readrawdata import readrawdata

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

#-------------------------------------------------------------------------------

def ll(filin,xm,xp):
    f = open(filin,"r")
    x = [] ; y = [] ; z = []
    lines = f.readlines()
    for iline in range(2,len(lines)):
        line = lines[iline].strip().split()
        xxx  = float(line[0])
        if (xxx >= xm and xp >= xxx):
            x.append(xxx)
            y.append(float(line[1])) 
            z.append(float(line[2])) 
    f.close()
    return x,y,z

#-------------------------------------------------------------------------------

def lp(filin,xm,xp):
    f = open(filin,"r")
    x = [] ; y = [] ; z = []
    lines = f.readlines()
    for iline in range(2,len(lines)):
        line = lines[iline].strip().split()
        if (len(line) < 6): break
        xxx  = float(line[3])
        if (xxx >= xm and xp >= xxx):
            x.append(xxx)
            y.append(float(line[4])) 
            z.append(float(line[5])) 
    f.close()
    return x,y,z

#-------------------------------------------------------------------------------

factor=27.211396132 ## conversion hartree -> eV

#-------------------------------------------------------------------------------
# Read arguments

xxmin = -1.e99
xxmax =  1.e99

narg = len(sys.argv)-1

if (narg > 0): infile = str(   sys.argv[1])
if (narg > 1): xxmin  = float( sys.argv[2])
if (narg > 2): xxmax  = float( sys.argv[3])
if (narg == 0 or narg > 3):
    print( "\n ERROR: Invalid number of arguments.\n")
    sys.exit(" Usage: PLOT-raman-epsilon.py FILENAME XMIN XMAX\n")
    
if ( not os.path.exists(infile)):
    sys.exit("\n ERROR: File "+infile+" does not exist!\n")

#-------------------------------------------------------------------------------
# Read data from infile

upho = [] ; epsr = [] ; epsi = []
upho, epsr, epsi = ll(infile,xxmin,xxmax)

uphop = [] ; epsrp = [] ; epsip = []
uphop, epsrp, epsip = lp(infile,xxmin,xxmax)

#-------------------------------------------------------------------------------
# Settings for the plot 
    
figcolor = 'white'

dpi = 300

fig = plt.figure(figsize=(18,10),dpi=dpi)

fig.patch.set_edgecolor(figcolor)
fig.patch.set_facecolor(figcolor)

mpl.rcParams['axes.linewidth']  = 4.0     # set the value globally
mpl.rcParams['grid.linewidth']  = 1.5
mpl.rcParams['xtick.labelsize'] = 30
mpl.rcParams['ytick.labelsize'] = 30
mpl.rcParams['axes.edgecolor']  = 'black'
mpl.rcParams['axes.labelsize']  = 40      # fontsize of the x any y labels
mpl.rcParams['axes.labelcolor'] = 'black'
mpl.rcParams['axes.axisbelow']  = 'True'  # whether axis gridlines and ticks are below
                                          # the axes elements (lines, text, etc)
mpl.rcParams['legend.fontsize'] = 30
plt.rcParams['xtick.major.pad'] = 10
plt.rcParams['ytick.major.pad'] = 10

#-------------------------------------------------------------------------------
# Plot 

ax1 = fig.add_subplot( 121)
ax2 = fig.add_subplot( 122, sharex=ax1)
 
#fig.subplots_adjust(hspace=-0.3)

ax1.set_xlabel('Phonon amplitude [Bohr]', labelpad=20)
ax2.set_xlabel('Phonon amplitude [Bohr]', labelpad=20)

ax1.set_ylabel(u'Re $\mathregular{\\varepsilon^{\\alpha\\beta}(\omega_L)}$', labelpad=20)
ax2.set_ylabel(u'Im $\mathregular{\\varepsilon^{\\alpha\\beta}(\omega_L)}$', labelpad=20)

for line in ax1.get_xticklines() + ax1.get_yticklines():
    line.set_markersize(10)
    line.set_markeredgewidth(3)

for line in ax2.get_xticklines() + ax2.get_yticklines():
    line.set_markersize(10)
    line.set_markeredgewidth(3)

ax1.plot( upho, epsr, color='red', lw=4.0, label="fit")
ax2.plot( upho, epsi, color='red', lw=4.0)

ax1.plot( uphop, epsrp, marker='o', color='blue', markersize=15.0, linewidth=0, label="data")
ax2.plot( uphop, epsip, marker='o', color='blue', markersize=15.0, linewidth=0)

leg=ax1.legend(loc=2,borderaxespad=0.7,numpoints=1)
leg.get_frame().set_linewidth(4.0)
leg.get_frame().set_edgecolor("grey")
leg.draw_frame(True)

if (narg > 1): plt.xlim( xmin=xxmin )
if (narg > 2): plt.xlim( xmax=xxmax )

ymin1, ymax1 = ax1.set_ylim()
ymin2, ymax2 = ax2.set_ylim()

ax1.vlines(0.0,ymin1,ymax1,linewidth=3.0,linestyles="dashed")
ax2.vlines(0.0,ymin2,ymax2,linewidth=3.0,linestyles="dashed")

fig.suptitle(infile,x=0.5,y=0.97,fontsize=40,color="darkblue")

fig.tight_layout(rect=[0, 0.01, 1, 0.92])
fig.savefig('PLOT.png',format='png',bbox_inches=0,dpi=300)
fig.savefig('PLOT.eps',format='eps',bbox_inches=0)

sys.exit()    

#-------------------------------------------------------------------------------
