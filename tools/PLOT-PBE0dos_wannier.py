#!/usr/bin/env python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

factor=27.211396132 ## conversion hartree -> eV


###################
# Read input data #
###################

## Create the list of input directories ##
root=os.getcwd()

# DOS from TDOS.OUT
dos0=[]
list1=[]
list2=[]
infile=root+"/TDOS.OUT"
for line in open(infile):
    i_line=line.split()
    if len(i_line):
       list1.append(float(i_line[0])*factor)
       list2.append(float(i_line[1])/factor)
    else:
       dos0.append([list1,list2])
       list1=[]
       list2=[]

# DOS from TDOS_WANNIER.OUT
dosw=[]
list1=[]
list2=[]
infile2=root+"/TDOS_WANNIER.OUT"
for line in open(infile2):
    i_line=line.split()
    if len(i_line):
       list1.append(float(i_line[0])*factor)
       list2.append(float(i_line[1])/factor) # convert to eV
    else:
       dosw.append([list1,list2])
       list1=[]
       list2=[]

################################################################################
## Settings for the plot #######################################################
################################################################################
    
figcolor = 'white'
dpi = 300
fig, ax1 = plt.subplots( 1, 1, figsize=(16,10), dpi=dpi)
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

#############################
##    Bandstructure plot   ##
#############################

ax1.xaxis.set_label_position('bottom')
ax1.yaxis.label.set_size(40)
ax1.set_ylabel('DOS [states/eV/unit cell]')
ax1.tick_params(axis='both', which='major', labelsize=mpl.rcParams['axes.labelsize'])


# Tick size
for line in ax1.get_xticklines() + ax1.get_yticklines():
    line.set_markersize(10)
    line.set_markeredgewidth(2)

ymax=max(max(dos0[0][1]),max(dosw[0][1]))
xmin=min(min(dos0[0][0]),min(dosw[0][0]))
xmax=max(max(dos0[0][0]),max(dosw[0][0]))
ax1.fill_between(dosw[0][0],dosw[0][1],color='r')
legend1 = mpatches.Patch( color='r', label='DOS Wannier')
legend2, = ax1.plot(dos0[0][0],dos0[0][1],'b',lw=2.0,label='DOS coarse grid')

leg=ax1.legend(handles=[legend1,legend2],bbox_to_anchor=(0.3,1),loc=2,borderaxespad=0.)
leg.draw_frame(True)

ax1.set_xlim(xmin,xmax)
ax1.set_ylim(0,1.2*ymax)
ax1.grid( True)
fig.savefig('PBE0_PBE_dos_wannier.png',format='png',bbox_inches=0,dpi=300)
fig.savefig('PBE0_PBE_dos_wannier.eps',format='eps',bbox_inches=0)

plt.show()
sys.exit()    
