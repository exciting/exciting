#!/usr/bin/env python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as pyl
import os

factor=27.211396132 ## conversion hartree -> eV


###################
# Read input data #
###################

## Create the list of input directories ##
root=os.getcwd()

print "\n################################################\n"
print "    Enter the working directories \n"
print "################################################\n"
LDA_dir=raw_input("name of PBE directory: ")
EXX_dir=raw_input("name of PBE0 directory: ")
print "\n"


# DFT states from DFT/BAND.OUT
ksene=[]
list1=[]
list2=[]
infile=root+"/"+LDA_dir+"/BAND.OUT"
for line in open(infile):
    i_line=line.split()
    if len(i_line):
       list1.append(float(i_line[0]))
       list2.append(float(i_line[1])*factor) # convert to eV
    else:
       ksene.append([list1,list2])
       list1=[]
       list2=[]

# EXX states from EXX/BAND.OUT
gwene=[]
list1=[]
list2=[]
infile2=root+"/"+EXX_dir+"/BAND.OUT"
for line in open(infile2):
    i_line=line.split()
    if len(i_line):
       list1.append(float(i_line[0]))
       list2.append(float(i_line[1])*factor) # convert to eV
    else:
       gwene.append([list1,list2])
       list1=[]
       list2=[]


# Read info about x-ticks position
bandlines=[]
infile3=root+"/"+EXX_dir+"/BANDLINES.OUT"
fid=open(infile3)
while 1:
    line=fid.readline()
    if not line:
        break
    i_line=line.split()
    bandlines.append(float(i_line[0]))
    # skip next two lines
    fid.readline()
    fid.readline()

# position of VBM (to be shifted to zero)
ivbm=4
ks0=max(ksene[ivbm-1][1])
gw0=max(gwene[ivbm-1][1])
for i in range(len(ksene)):
    for j in range(len(ksene[i][1])):
        ksene[i][1][j]=ksene[i][1][j]-ks0

for i in range(len(gwene)):
    for j in range(len(gwene[i][1])):
        gwene[i][1][j]=gwene[i][1][j]-gw0

################################################################################
## Settings for the plot #######################################################
################################################################################
    
figcolor = 'white'
dpi = 100
fig = plt.figure(figsize=(16,10),dpi=dpi)
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

ax1 = fig.add_axes([0.1,0.1,0.8,0.8])
ax1.xaxis.grid(True,which='major',color='k',linestyle='-',linewidth=2)
ax1.xaxis.set_label_position('bottom')
ax1.set_xticks(bandlines)
labels = ax1.set_xticklabels(('W','L','$\Gamma$','X','W','K'))
ax1.set_ylabel('Energy [eV]')


# Tick size
for line in ax1.get_xticklines() + ax1.get_yticklines():
    line.set_markersize(10)
    line.set_markeredgewidth(2)

ymin=1000.0
ymax=-1000.0
## get values for ymin and ymax
len1=len(ksene) 
len2=len(gwene)
if len1<len2:
        bandlen=len1
	for i in range(len(ksene)):
		y=min(ksene[i][1])
    		if (y<ymin):
        		ymin=y
    		y=max(ksene[i][1])
    		if (y>ymax):
        		ymax=y

else:
        bandlen=len2
        for i in range(len(gwene)):
                y=min(gwene[i][1])
                if (y<ymin):
                        ymin=y
                y=max(gwene[i][1])
                if (y>ymax):
                        ymax=y


for i in range(bandlen-1):
    ax1.plot(ksene[i][0],ksene[i][1],'b',lw=3.0)
    ax1.plot(gwene[i][0],gwene[i][1],'r',lw=3.0)
i=bandlen-1
ax1.plot(ksene[i][0],ksene[i][1],'b',lw=3.0,label='PBE')
ax1.plot(gwene[i][0],gwene[i][1],'r',lw=3.0,label='PBE0')

leg=ax1.legend(bbox_to_anchor=(0.825,0.23),loc=2,borderaxespad=0.)
leg.draw_frame(True)

# add zero level
x0=[0.0,max(ksene[0][0])]
y0=[0.0,0.0]
ax1.plot(x0,y0,'k:',lw=1.0)
#ymax=40.0
ax1.set_xlim(0,max(ksene[0][0]))
ax1.set_ylim(-25,15)
pyl.grid(True)
fig.savefig('PBE0_PBE.png',format='png',bbox_inches=0,dpi=45)
fig.savefig('PBE0_PBE.eps',format='eps',bbox_inches=0)

plt.show()
sys.exit()    
