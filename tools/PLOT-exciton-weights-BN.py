#!/usr/bin/env python
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

ha2ev = 27.211396132 ## conversion hartree -> eV
bo2an = 0.52917721067

narg = len(sys.argv)-1
if( narg < 1):
    print( "\n** ERROR: Must specify name of file and direction on command line.\n")
    print( "** Usage: ", sys.argv[0], " <KS band structure width weights>")
    sys.exit(0)
if( not os.path.isfile( sys.argv[1])):
    print( "\n** ERROR: Input file %s was not found." % sys.argv[1])
    sys.exit(0)
inname = sys.argv[1].strip()
if( narg == 1):
    scale = 1.0
elif( narg == 2):
    scale = float( sys.argv[2])
elif( narg == 3): 
    ymin = float( sys.argv[2])
    ymax = float( sys.argv[3])
    scale = 1.0
else:
    ymin = float( sys.argv[2])
    ymax = float( sys.argv[3])
    scale = float( sys.argv[4])
scale *= 2e2
tmp = inname.split( '.')
tmp[-1] = 'png'
outname1 = '.'.join( tmp)
tmp[-1] = 'pdf'
outname2 = '.'.join( tmp)

###################
# Read input data #
###################
bandwgt, bandwgtdim = readrawdata( inname)
bandlin, bandlindim = readrawdata( '../BANDLINES.OUT')
bandlines = bandlin[:,0,0]

#########################
# Settings for the plot #
#########################
    
figcolor = 'white'
fig = plt.figure( figsize=(10,10))
fig.patch.set_edgecolor(figcolor)
fig.patch.set_facecolor(figcolor)

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

######################
# Bandstructure plot #
######################

ax1 = fig.add_axes( [0.17,0.1,0.75,0.8])
ax2 = ax1.twinx()
ax1.xaxis.grid( True, which='major', color='k', linestyle='-', linewidth=2)
ax1.xaxis.set_label_position( 'bottom')
ax1.set_xticks( bandlines)
ax1.set_xticklabels( ( 'X', 'W','L','$\Gamma$','X'))
ax1.set_ylabel( 'Energy [eV]')

# Tick size
for line in ax1.get_xticklines() + ax1.get_yticklines():
    line.set_markersize(10)
    line.set_markeredgewidth(2)

xmin = np.amin( bandwgt[:,0,:])
xmax = np.amax( bandwgt[:,0,:])
ax1.set_xlim( [xmin, xmax])
if( narg > 2): ax1.set_ylim( [ymin, ymax])

# band-structure and weights
for i in range( bandwgtdim[2]):
    ax1.plot( bandwgt[i,0,:], bandwgt[i,1,:], 'b', lw=3.0, zorder=10)
for i in range( bandwgtdim[2]):
    ax1.scatter( bandwgt[i,0,:], bandwgt[i,1,:], s=(scale*bandwgt[i,2,:])**2, lw=3.0, edgecolor='r', facecolor='none', zorder=11)

# Fermi level
ax1.plot( [xmin, xmax], [0, 0], 'k', lw=3.0, ls='-')
ax2.set_ylim( ax1.get_ylim())
ax2.set_yticks( [0, 0])
ax2.set_yticklabels( ('$\\mathregular{E_{F}}$', ''))

ax1.grid( True)
plt.title( "BN excitonic weights", fontsize=mpl.rcParams['ytick.labelsize'], y=1.03)

fig.savefig( outname1, format='png', bbox_inches=0, dpi=300)
fig.savefig( outname2, format='pdf', bbox_inches=0)

plt.show()
sys.exit()    
