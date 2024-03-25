#!/usr/bin/python
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
import re

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

#-------------------------------------------------------------------------------
# Check arguments

yylabel = r"Im $\mathregular{\epsilon_M}$"
narg = len(sys.argv)-1
nfiles = narg
function = 2
n0 = 0

if narg<1:
    print "\nERROR: Nothing to plot!\n"
    print "**Usage**:    PLOT-spectra.py [--imag/--real] file1.xml [ file2.xml [ file3.xml [...]]]\n"
    sys.exit()

if ( str(sys.argv[1]).lower() == "--real" ):
    function = 1
    nfiles = narg - 1
    n0 = 1
    yylabel = r"Re $\mathregular{\epsilon_M}$"
if ( str(sys.argv[1]).lower() == "--imag" ):
    function = 2
    nfiles = narg - 1
    n0 = 1
    yylabel = r"Im $\mathregular{\epsilon_M}$"

fnames=[]
for i in range(n0,narg):
    fnames.append(sys.argv[i+1])
    if not os.path.isfile(fnames[i-n0]):
        print "Error: file \"%s\" doesn't exist"%(fnames[i])
        sys.exit()

#-------------------------------------------------------------------------------
# Parse EPSILON/LOSS function data files

xdata=[] ; ydata=[] ; labels=[] ; legends=[]
print

for i,fname in enumerate(fnames):
    xdata.append([])
    ydata.append([])
    labels.append({})
    print " Parsing "+fname
    sfname = fname.split("/")[-1]
    sfname = re.split('_|-',sfname)

    if "BSE"          in sfname: legend="BSE "
    if "TDA"          in sfname: legend=legend+"TDA "
    if  not "TDA"     in sfname and "BSE" in sfname: legend=legend+"full "
    if "singlet"      in sfname: legend=legend+"singlet"
    if "triplet"      in sfname: legend=legend+"triplet"
    if "RPA"          in sfname: legend=legend+"RPA"
    if "IP"           in sfname: legend=legend+"IP"
    if "FXCRPA"       in sfname: legend="RPA "
    if "FXCALDA"      in sfname: legend="ALDA "
    if "FXCLRCstatic" in sfname: legend="LRCstatic "
    if "FXCRBO"       in sfname: legend="RBO "
    if "FXCLRCdyn"    in sfname: legend="LRCdyn "
    if "FXCMB1"       in sfname: legend="MB1 "
    if  not "NLF"     in sfname and not "BSE" in sfname: legend=legend+"(LFE) "
    if "NLF"          in sfname: legend=legend+"(no-LFE) "
    if sfname[-2][0:2]=="OC" and "LOSS" in sfname: legend=legend+"Optical(%s)"%(sfname[-2][2:])

    legends.append(legend)

    if "LOSS" in sfname:
        labels[i]["ylabel"] = "Loss function"
    if "EPSILON" in sfname:
        labels[i]["ylabel"] = "Im $\epsilon_M$"

    labels[i]["xlabel"] = "Energy"

    file = open(fname,'r')
    for lines in file:
        if 'Frequency' in lines:
            break

    for lines in file:
        data = lines.split()
        xdata[i].append(float(data[0]))
        ydata[i].append(float(data[function]))
    file.close()

print

#-------------------------------------------------------------------------------
# Settings for the plot

figcolor = 'white'
dpi = 300
fig = plt.figure( figsize=(15,10), dpi=dpi)
fig.patch.set_edgecolor(figcolor)
fig.patch.set_facecolor(figcolor)

mpl.rcParams['axes.linewidth'] = 4.0 # set the value globally
mpl.rcParams['grid.linewidth'] = 1.5
mpl.rcParams['xtick.labelsize'] = 30
mpl.rcParams['ytick.labelsize'] = 30
mpl.rcParams['axes.edgecolor'] = 'black'
mpl.rcParams['axes.labelsize'] = '40'     # fontsize of the x any y labels
mpl.rcParams['axes.labelcolor'] = 'black'
mpl.rcParams['axes.axisbelow'] = 'True'   # whether axis gridlines and ticks are below
                                          # the axes elements (lines, text, etc)
mpl.rcParams['legend.fontsize'] = '40'
plt.rcParams['xtick.major.pad'] = '10'
plt.rcParams['ytick.major.pad'] = '10'

colors=['firebrick','mediumblue','g','y','c','m','k']

#-------------------------------------------------------------------------------
# Plot

ax1 = fig.add_subplot( 111)
ax1.xaxis.set_label_position('bottom')
ax1.set_xlabel( 'Energy [eV]')
if "EPSILON" in sfname: ax1.axhline( y=0, linestyle="dashed", linewidth=3, color="black")

# Tick size
for line in ax1.get_xticklines() + ax1.get_yticklines():
    line.set_markersize(10)
    line.set_markeredgewidth(2)

hndls = []
xmin = 1e100 ; xmax = -1e100
ymin = 1e100 ; ymax = -1e100
for i in range(nfiles):
    h, = ax1.plot( xdata[i], ydata[i], 'r-', color=colors[np.mod(i,7)], label=legends[i], lw=3.0)
    hndls.append( h)
    xmin = min( xmin, min(xdata[i]) )
    xmax = max( xmax, max(xdata[i]) )
    ymin = min( ymin, min(ydata[i]) )
    ymax = max( ymax, max(ydata[i]) )

leg = ax1.legend( handles=hndls, loc='best', borderaxespad=0.5)
leg.get_frame().set_edgecolor('grey')
leg.get_frame().set_linewidth( 4.0)
leg.draw_frame(True)

ax1.set_xlim( 0.0, xmax)
if "EPSILON" in sfname: ax1.set_ylim( ymin-0.05*(ymax-ymin), ymax+0.05*(ymax-ymin) )
if "LOSS"    in sfname: ax1.set_ylim( ymin, ymax+0.2*(ymax-ymin) )
ax1.set_xlabel(str.capitalize(labels[0]["xlabel"])+" [eV]")
if "EPSILON" in sfname:
    ax1.set_ylabel(yylabel, labelpad=20)
else:
    ax1.set_ylabel(str.capitalize(labels[0]["ylabel"]).encode('string-escape'), labelpad=20)

fig.tight_layout()
plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=300)

#-------------------------------------------------------------------------------
