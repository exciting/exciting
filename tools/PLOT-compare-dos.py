#!/usr/bin/env python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   lxml  import etree as et
import sys
import matplotlib         as mpl
import matplotlib.pyplot  as plt
import matplotlib.patches as mpatches
import pylab              as pyl
import os
import numpy              as np
from readrawdata import readrawdata

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

factor=27.211396132 ## conversion hartree -> eV

#-------------------------------------------------------------------------------
# Read arguments

narg = len(sys.argv)-1
lgw  = False
lwa  = False
case = " "

if( narg == 1):
    case = str( sys.argv[1])
elif( narg == 2):
    xxmin = float( sys.argv[1])
    xxmax = float( sys.argv[2])
elif( narg == 3):
    xxmin = float( sys.argv[1])
    xxmax = float( sys.argv[2])
    case = str( sys.argv[3])
elif( narg > 3):
    print( "\n ERROR: Invalid number of arguments.\n")
    sys.exit(" Usage: PLOT-compare-dos.py Energy_min Energy_max CASE\n")

if (case.strip().lower() == "gw"): lgw = True
if (case.strip().lower() == "wannier"): lwa = True

if lgw: print( "\n Comparison for GW calculations\n")
if lwa: print( "\n Comparison for WANNIER calculations\n")

if ( os.path.exists('GW_INFO.OUT') and not lgw ):
    print( "\n ERROR: This is a GW directory! Check your command line:\n")
    sys.exit(" Usage: PLOT-compare-dos.py Energy_min Energy_max GW\n")
    
if ( not os.path.exists('GW_INFO.OUT') and lgw ):
    sys.exit(" ERROR: This is NOT a GW directory! Delete GW from the command line!\n")

if ( os.path.exists('WANNIER_INFO.OUT') and not lwa ):
    print( "\n ERROR: This is a WANNIER directory! Check your command line:\n")
    sys.exit(" Usage: PLOT-compare-dos.py Energy_min Energy_max WANNIER\n")
    
if ( not os.path.exists('WANNIER_INFO.OUT') and lwa ):
    sys.exit(" ERROR: This is NOT a WANNIER directory! Delete WANNIER from the command line!\n")

#-------------------------------------------------------------------------------
# Create the list of input directories 

root=os.getcwd()

if(not lgw and not lwa):
    print( "\n################################################\n")
    print( " Enter the names of the 2 directories to compare\n")
    print( "------------------------------------------------\n")
    dirone=raw_input(" Directory 1 ==> ")
    dirtwo=raw_input(" Directory 2 ==> ")
    print 
    print( "################################################\n")

#-------------------------------------------------------------------------------
# Read data from dirone

if lgw:
    infile=root+"/TDOS.OUT"
    label1 = "KS"
elif lwa:
    infile=root+"/TDOS.OUT"
    label1 = "original grid"
else:
    infile=root+"/"+dirone+"/TDOS.OUT"
    label1 = dirone
dos1, dim1 = readrawdata( infile)
dos1[:,0,:] *= factor
dos1[:,1,:] /= factor

#-------------------------------------------------------------------------------
# Read data from dirtwo

if lgw:
    infile=root+"/TDOS-QP.OUT"
    label2 = "$\mathregular{G_0W_0}$"
elif lwa:
    infile=root+"/TDOS_WANNIER.OUT"
    label2 = "Wannier"
else:
    infile = root+"/"+dirtwo+"/TDOS.OUT"
    label2 = dirtwo
dos2, dim2 = readrawdata( infile)
dos2[:,0,:] *= factor
dos2[:,1,:] /= factor

#-------------------------------------------------------------------------------
# Settings for the plot 
    
figcolor = 'white'
dpi = 300
fig = plt.figure(figsize=(15,10),dpi=dpi)
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
# DOS plot 

#ax1 = fig.add_axes([0.14,0.1,0.8,0.8])
ax1 = fig.add_subplot( 111)

ax1.xaxis.set_label_position('bottom')
ax1.set_xlabel('Energy [eV]')
ax1.set_ylabel('DOS [states/eV/unit cell]', labelpad=20)

for line in ax1.get_xticklines() + ax1.get_yticklines():
    line.set_markersize(10)
    line.set_markeredgewidth(2)

ax1.fill_between( dos2[0,0,:], dos2[0,1,:], color='darksalmon')
ax1.plot( dos2[0,0,:], dos2[0,1,:], color='firebrick', lw=3.0)
legend2 = mpatches.Patch( facecolor='darksalmon', edgecolor='firebrick', lw=3.0, label=label2)
legend1, = ax1.plot( dos1[0,0,:], dos1[0,1,:], color='mediumblue', lw=3.0, label=label1)

leg = ax1.legend( handles=[legend1, legend2], loc='best',borderaxespad=0.5)
leg.get_frame().set_linewidth(4.0)
leg.get_frame().set_edgecolor("grey")
leg.draw_frame(True)

if ( narg < 2):
    xxmin = max( np.amin( dos1[0,0,:]), np.amin( dos2[0,0,:]))
    xxmax = min( np.amax( dos1[0,0,:]), np.amax( dos2[0,0,:]))
plt.xlim( xmin=xxmin, xmax=xxmax)

fig.tight_layout()
fig.savefig('PLOT.png',format='png',bbox_inches=0,dpi=300)
fig.savefig('PLOT.eps',format='eps',bbox_inches=0)

sys.exit()    

#-------------------------------------------------------------------------------
