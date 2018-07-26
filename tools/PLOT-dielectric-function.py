#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

import sys
import matplotlib         as mpl
import matplotlib.pyplot  as plt
import numpy              as np
import os
from readrawdata import readrawdata

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

#-------------------------------------------------------------------------------
# Read arguments

narg = len(sys.argv)-1
opt = 'i'
if( narg < 1):
    print( "\n** ERROR: Must specify name of file and direction on command line.\n")
    print( "** Usage: ", sys.argv[0], " <dielectric function>")
    sys.exit(0)
inname = sys.argv[1].strip()
if( not os.path.isfile( inname)):
    print( "\n** ERROR: Input file %s was not found." % sys.argv[1])
    sys.exit(0)
if( narg > 1):
    if( sys.argv[2].strip() == 'i'):
        opt = 'i'
    elif( sys.argv[2].strip() == 'r'):
        opt = 'r'
    elif( sys.argv[2].strip() == 'b'):
        opt = 'b'
    else:
        print( "\n** ERROR: Input parameter 2 is invalid. Use either 'r' (real part), 'i' (imaginary part, default) or 'b' (both).")
        sys.exit(0)

#-------------------------------------------------------------------------------
# Read data

diel, dieldim = readrawdata( inname)

#-------------------------------------------------------------------------------
# Settings for the plot

figcolor = 'white'
fig = plt.figure( figsize=(12,8))
fig.figurePatch.set_edgecolor(figcolor)
fig.figurePatch.set_facecolor(figcolor)

mpl.rcParams['axes.linewidth'] = 3.0 # set the value globally
mpl.rcParams['grid.linewidth'] = 1.5
mpl.rcParams['xtick.labelsize'] = 30
mpl.rcParams['ytick.labelsize'] = 30
mpl.rcParams['axes.edgecolor'] = 'black'
mpl.rcParams['axes.labelsize'] = '30'     # fontsize of the x any y labels
mpl.rcParams['axes.labelcolor'] = 'black'
mpl.rcParams['axes.axisbelow'] = 'True'   # whether axis gridlines and ticks are below
                                          # the axes elements (lines, text, etc)
mpl.rcParams['legend.fontsize'] = '25'
plt.rcParams['xtick.major.pad'] = '10'
plt.rcParams['ytick.major.pad'] = '10'

#-------------------------------------------------------------------------------
# Plot

ax1 = fig.add_axes( [0.15,0.15,0.8,0.75])
ax1.set_xlabel( 'Energy [eV]')
tmp = inname.split( '.')
if( opt == 'r'):
    ax1.set_ylabel( 'Re $\\mathregular{\\epsilon_{M}}$', labelpad=20)
    tmp[-2] += '_real'
elif( opt == 'i'):
    ax1.set_ylabel( 'Im $\\mathregular{\\epsilon_{M}}$', labelpad=20)
    tmp[-2] += '_imag'
elif( opt == 'b'):
    tmp[-2] += '_full'
    ax1.set_ylabel( '$\\mathregular{\\epsilon_{M}}$', labelpad=20)
tmp[-1] = 'png'
outname1 = '.'.join( tmp)
tmp[-1] = 'pdf'
outname2 = '.'.join( tmp)

# Tick size
for line in ax1.get_xticklines() + ax1.get_yticklines():
    line.set_markersize(10)
    line.set_markeredgewidth(2)

xmin = np.amin( diel[:,0,:])
xmax = np.amax( diel[:,0,:])
ax1.set_xlim( [xmin, xmax])

if( opt == 'r'):
    ax1.plot( diel[0,0,:], diel[0,1,:], 'r-', label="real", color="royalblue", lw=3.0)
elif( opt == 'i'):
    ax1.plot( diel[0,0,:], diel[0,2,:], 'r-', label="imaginary", color="firebrick", lw=3.0)
elif( opt == 'b'):
    leg1, = ax1.plot( diel[0,0,:], diel[0,1,:], 'r-', label="real", color="royalblue", lw=3.0)
    leg2, = ax1.plot( diel[0,0,:], diel[0,2,:], 'r-', label="imaginary", color="firebrick", lw=3.0)
    leg = ax1.legend( handles=[leg1,leg2], loc='best', borderaxespad=1.0)
    leg.get_frame().set_edgecolor('#CCCCCC')
    leg.get_frame().set_linewidth( 3)
		 
ax1.grid( True)
plt.title( "Macroscopic dielectric function", fontsize=mpl.rcParams['ytick.labelsize'], y=1.03)

fig.savefig( outname1, format='png', bbox_inches='tight', dpi=300)
fig.savefig( outname2, format='pdf', bbox_inches=0)

plt.show()
sys.exit()    

#-------------------------------------------------------------------------------
