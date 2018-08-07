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
    yymin = float( sys.argv[1])
    yymax = float( sys.argv[2])
elif( narg == 3):
    yymin = float( sys.argv[1])
    yymax = float( sys.argv[2])
    case = str( sys.argv[3])
elif( narg == 0 or narg > 3):
    print( "\n ERROR: Invalid number of arguments.\n")
    sys.exit(" Usage: PLOT-compare-bands.py Energy_min Energy_max CASE\n")

if (case.strip().lower() == "gw"): lgw = True
if (case.strip().lower() == "wannier"): lwa = True

if lgw: print( "\n Comparison for GW calculations\n")
if lwa: print( "\n Comparison for WANNIER calculations\n")

if ( os.path.exists('GW_INFO.OUT') and not lgw ):
    print( "\n ERROR: This is a GW directory! Check your command line:\n")
    sys.exit(" Usage: PLOT-compare-bands.py Energy_min Energy_max GW\n")
    
if ( not os.path.exists('GW_INFO.OUT') and lgw ):
    sys.exit(" ERROR: This is NOT a GW directory! Delete GW from the command line!\n")

#-------------------------------------------------------------------------------
# Create the list of input directories 

root=os.getcwd()

if (not lgw):
    print( "\n################################################\n")
    print( " Enter the names of the 2 directories to compare\n")
    print( "------------------------------------------------\n")
    dirone=raw_input(" Directory 1 ==> ")
    dirtwo=raw_input(" Directory 2 ==> ")
    print 
    print( "################################################\n")

    wandir = 0
    if ( os.path.exists(root+"/"+dirone+"/WANNIER_INFO.OUT")):
        wandir = 1
        if( not lwa):
            print( "\n ERROR: Directory 1 is a WANNIER directory! Check your command line:\n")
            sys.exit(" Usage: PLOT-compare-bands.py Energy_min Energy_max WANNIER\n")
        
    if ( os.path.exists(root+"/"+dirtwo+"/WANNIER_INFO.OUT")):
        if( wandir == 1):
            sys.exit( "\n ERROR: Both directories are WANNIER directories! Nothing to compare with.\n")
        else:
            wandir = 2
        if( not lwa):
            print( "\n ERROR: Directory 2 is a WANNIER directory! Check your command line:\n")
            sys.exit(" Usage: PLOT-compare-bands.py Energy_min Energy_max WANNIER\n")
    
    if ( wandir == 0 and lwa):
        sys.exit(" ERROR: No directory is a WANNIER directory! Delete WANNIER from the command line!\n")

#-------------------------------------------------------------------------------
# Read data from dirone

if lgw:
    infile=root+"/BAND.OUT"
    label1 = "KS"
elif lwa:
    infile=root+"/"+dirone+"/BAND.OUT"
    label1 = dirone
    if( wandir == 1):
        label1 += " F"
else:
    infile=root+"/"+dirone+"/BAND.OUT"
    label1 = dirone
band1, dim1 = readrawdata( infile)
band1[:,1,:] *= factor

#-------------------------------------------------------------------------------
# Read data from dirtwo

if lgw:
    infile=root+"/BAND-QP.OUT"
    label2 = "$\mathregular{G_0W_0}$"
elif lwa:
    infile=root+"/"+dirtwo+"/BAND.OUT"
    label2 = dirtwo
    if( wandir == 2):
        label2 += " F"
else:
    infile = root+"/"+dirtwo+"/BAND.OUT"
    label2 = dirtwo
band2, dim2 = readrawdata( infile)
band2[:,1,:] *= factor

#-------------------------------------------------------------------------------
# Read Wannier bands if needed

if lwa:
    if( wandir == 1):
        infile=root+"/"+dirone+"/BAND_WANNIER.OUT"
        label3 = dirone+" W"
    else:
        infile=root+"/"+dirtwo+"/BAND_WANNIER.OUT"
        label3 = dirtwo+" W"
    band3, dim3 = readrawdata( infile)
    band3[:,1,:] *= factor

#-------------------------------------------------------------------------------
# Read info about x-ticks position

if lgw:
    infile=root+"/BANDLINES.OUT"
else:
    infile=root+"/"+dirtwo+"/BANDLINES.OUT"
bandlin, bandlindim = readrawdata( infile)
bandlines = bandlin[:,0,0]
    
#-------------------------------------------------------------------------------
# Read info about x-ticks labels
  
if lgw:
    infile = root+"/input.xml"
else:
    infile = root+"/"+dirtwo+"/input.xml"
ifile = open( infile, "r")
iroot = et.parse(ifile).getroot()
iphod = -1
while True: 
    iphod = iphod+1
    if ( str(iroot[iphod].tag) == "properties" ): break

iphpl = -1
while True: 
    iphpl = iphpl+1
    if ( str(iroot[iphod][iphpl].tag) == "bandstructure" ): break

ipath = iroot[iphod][iphpl][0][0]
llist = []
for i in range(len(ipath)):
    label = " "
    if 'label' in ipath[i].attrib: label = ipath[i].attrib['label']
    if ( label.strip().lower() == 'gamma'): label = u'\u0393'
    llist.append(label)

#-------------------------------------------------------------------------------
# If not GW, set zero to the value of the VBM (for both directories)

if (not lgw):
    ivbm=4
    bnd0 = np.amax( band1[ivbm-1,1,:])
    band1[:,1,:] -= bnd0
    bnd0 = np.amax( band2[ivbm-1,1,:])
    band2[:,1,:] -= bnd0
    if lwa: 
        bnd0 = np.amax( band3[ivbm-1,1,:])
        band3[:,1,:] -= bnd0
#-------------------------------------------------------------------------------
# Settings for the plot 
    
figcolor = 'white'
dpi = 300
if lwa:
    fig = plt.figure(figsize=(25,10),dpi=dpi)
else:
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
# Band structure plot 

#ax1 = fig.add_axes([0.14,0.1,0.8,0.8])
if lwa:
    ax1 = fig.add_subplot( 121)
    ax2 = fig.add_subplot( 122, sharex=ax1, sharey=ax1)
    ax1.xaxis.grid(True,which='major',color='k',linestyle='-',linewidth=3)
    ax2.xaxis.grid(True,which='major',color='k',linestyle='-',linewidth=3)
    ax1.axhline(y=0,linestyle="dashed",linewidth=3,color="black")
    ax2.axhline(y=0,linestyle="dashed",linewidth=3,color="black")
    plt.setp( ax2.get_yticklabels(), visible=False)
else:
    ax1 = fig.add_subplot( 111)
    ax1.xaxis.grid(True,which='major',color='k',linestyle='-',linewidth=3)
    ax1.axhline(y=0,linestyle="dashed",linewidth=3,color="black")

ax1.xaxis.set_label_position('bottom')
ax1.set_xticks(bandlines)
ax1.set_xticklabels(llist)
ax1.set_ylabel('Energy [eV]', labelpad=20)

for line in ax1.get_xticklines() + ax1.get_yticklines():
    line.set_markersize(10)
    line.set_markeredgewidth(2)
if lwa:
    for line in ax2.get_xticklines() + ax2.get_yticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(2)

bandlen = min( dim1[2], dim2[2])
if lwa:
    bandlen = min( bandlen, dim3[2])
                        
for i in range(bandlen-1):
    ax1.plot( band1[i,0,:], band1[i,1,:], color='mediumblue', lw=3.0)
    ax1.plot( band2[i,0,:], band2[i,1,:], color='firebrick', lw=3.0)
i=bandlen-1

ax1.plot( band1[i,0,:], band1[i,1,:], color='mediumblue', lw=3.0, label=label1)
ax1.plot( band2[i,0,:], band2[i,1,:], color='firebrick', lw=3.0, label=label2)

leg=ax1.legend(loc=4,borderaxespad=0.5)
leg.get_frame().set_linewidth(4.0)
leg.get_frame().set_edgecolor("grey")
leg.draw_frame(True)

if lwa:
    for i in range(bandlen-1):
        ax2.plot( band1[i,0,:], band1[i,1,:], color='mediumblue', lw=3.0)
        ax2.plot( band3[i,0,:], band3[i,1,:], color='firebrick', lw=3.0)
    i=bandlen-1
    
    ax2.plot( band1[i,0,:], band1[i,1,:], color='mediumblue', lw=3.0, label=label1)
    ax2.plot( band3[i,0,:], band3[i,1,:], color='firebrick', lw=3.0, label=label3)

    leg=ax2.legend(loc=4,borderaxespad=0.5)
    leg.get_frame().set_linewidth(4.0)
    leg.get_frame().set_edgecolor("grey")
    leg.draw_frame(True)

if ( narg > 1): plt.ylim( ymin=yymin, ymax=yymax)

fig.tight_layout()
fig.savefig('PLOT.png',format='png',bbox_inches=0,dpi=300)
fig.savefig('PLOT.eps',format='eps',bbox_inches=0)

sys.exit()    

#-------------------------------------------------------------------------------
