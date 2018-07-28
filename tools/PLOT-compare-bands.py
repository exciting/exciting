#!/usr/bin/env python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   lxml  import etree as et
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as pyl
import os

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

factor=27.211396132 ## conversion hartree -> eV

#-------------------------------------------------------------------------------
# Read arguments

narg = len(sys.argv)-1
lgw  = False
gw   = " "

if (narg > 0): yymin = float(sys.argv[1])
if (narg > 1): yymax = float(sys.argv[2])
if (narg > 2): gw    =   str(sys.argv[3])
if (gw == "GW" or gw == "gw"): lgw = True

if lgw: print "\n Comparison for GW calculations\n"

if ( os.path.exists('GW_INFO.OUT') and not lgw ):
    print "\n ERROR: This is a GW directory! Check your command line:\n"
    sys.exit(" Usage: PLOT-compare-bands.py Energy_min Energy_max GW\n")
    
if ( not os.path.exists('GW_INFO.OUT') and lgw ):
    sys.exit(" ERROR: This is NOT a GW directory! Delete GW from the command line!\n")

#-------------------------------------------------------------------------------
# Read input data

#-------------------------------------------------------------------------------
# Create the list of input directories 

root=os.getcwd()

if (not lgw):
    print "\n################################################\n"
    print " Enter the names of the 2 directories to compare\n"
    print "------------------------------------------------\n"
    dirone=raw_input(" Directory 1 ==> ")
    dirtwo=raw_input(" Directory 2 ==> ")
    print 
    print "################################################\n"

#-------------------------------------------------------------------------------
# Read data from dirone

ksene=[]
list1=[]
list2=[]
if  lgw:
    infile="BAND.OUT"
    dirone="KS"
else:
    infile=root+"/"+dirone+"/BAND.OUT"
for line in open(infile):
    i_line=line.split()
    if len(i_line):
       list1.append(float(i_line[0]))
       list2.append(float(i_line[1])*factor) # convert to eV
    else:
       ksene.append([list1,list2])
       list1=[]
       list2=[]

#-------------------------------------------------------------------------------
# Read data from dirtwo

gwene=[]
list1=[]
list2=[]
if  lgw:
    infile2="BAND-QP.OUT"
    dirtwo="$\mathregular{G_0W_0}$"
else:
    infile2=root+"/"+dirtwo+"/BAND.OUT"
for line in open(infile2):
    i_line=line.split()
    if len(i_line):
       list1.append(float(i_line[0]))
       list2.append(float(i_line[1])*factor) # convert to eV
    else:
       gwene.append([list1,list2])
       list1=[]
       list2=[]

#-------------------------------------------------------------------------------
# Read info about x-ticks position

bandlines=[]

if  lgw:
    infile3="BANDLINES.OUT"
else:
    infile3=root+"/"+dirtwo+"/BANDLINES.OUT"
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
    
#-------------------------------------------------------------------------------
# Read info about x-ticks labels
  
if  lgw:
    ifile = open("input.xml","r")
else:
    ifile = open(dirtwo+"/input.xml","r")
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
    if ( label == 'Gamma'): label = u'\u0393'
    if ( label == 'gamma'): label = u'\u0393'
    if ( label == 'GAMMA'): label = u'\u0393'
    llist.append(label)

#-------------------------------------------------------------------------------
# If not GW, set zero to the value of the VBM (for both directories)

if (not lgw):
    ivbm=4
    ks0=max(ksene[ivbm-1][1])
    gw0=max(gwene[ivbm-1][1])
    for i in range(len(ksene)):
        for j in range(len(ksene[i][1])):
            ksene[i][1][j]=ksene[i][1][j]-ks0

    for i in range(len(gwene)):
        for j in range(len(gwene[i][1])):
            gwene[i][1][j]=gwene[i][1][j]-gw0
        
#-------------------------------------------------------------------------------
# Settings for the plot 
    
figcolor = 'white'
dpi = 300
fig = plt.figure(figsize=(15,10),dpi=dpi)
fig.figurePatch.set_edgecolor(figcolor)
fig.figurePatch.set_facecolor(figcolor)

mpl.rcParams['axes.linewidth']  = 4.0     # set the value globally
mpl.rcParams['grid.linewidth']  = 1.5
mpl.rcParams['xtick.labelsize'] = 30
mpl.rcParams['ytick.labelsize'] = 30
mpl.rcParams['axes.edgecolor']  = 'black'
mpl.rcParams['axes.labelsize']  = 50      # fontsize of the x any y labels
mpl.rcParams['axes.labelcolor'] = 'black'
mpl.rcParams['axes.axisbelow']  = 'True'  # whether axis gridlines and ticks are below
                                          # the axes elements (lines, text, etc)
mpl.rcParams['legend.fontsize'] = 30
plt.rcParams['xtick.major.pad'] = 10
plt.rcParams['ytick.major.pad'] = 10

#-------------------------------------------------------------------------------
# Band structure plot 

ax1 = fig.add_axes([0.14,0.1,0.8,0.8])
ax1.xaxis.grid(True,which='major',color='k',linestyle='-',linewidth=3)
ax1.xaxis.set_label_position('bottom')
ax1.set_xticks(bandlines)
ax1.set_xticklabels(llist)
ax1.set_ylabel('Energy [eV]')

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
                        
ax1.axhline(y=0,linestyle="dashed",linewidth=3,color="black")

for i in range(bandlen-1):
    ax1.plot(ksene[i][0],ksene[i][1],'b',lw=3.0)
    ax1.plot(gwene[i][0],gwene[i][1],'r',lw=3.0)
i=bandlen-1

ax1.plot(ksene[i][0],ksene[i][1],'b',lw=3.0,label=dirone)
ax1.plot(gwene[i][0],gwene[i][1],'r',lw=3.0,label=dirtwo)

leg=ax1.legend(loc=4,borderaxespad=0.5)
leg.get_frame().set_linewidth(4.0)
leg.get_frame().set_edgecolor("grey")
leg.draw_frame(True)

#ax1.set_xlim(0,max(ksene[0][0]))
#ax1.set_ylim(-25,15)

if ( narg > 0): plt.ylim(ymin=yymin)
if ( narg > 1): plt.ylim(ymax=yymax)

fig.savefig('PLOT.png',format='png',bbox_inches=0,dpi=300)
fig.savefig('PLOT.eps',format='eps',bbox_inches=0)

sys.exit()    

#-------------------------------------------------------------------------------
