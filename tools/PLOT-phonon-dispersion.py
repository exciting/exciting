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

if (narg > 0): yymin  = float(sys.argv[1])
if (narg > 1): yymax  = float(sys.argv[2])
if (narg > 2): 
    plabel = str(sys.argv[3])
    mcolor = "b" 
   
#-------------------------------------------------------------------------------
# Read data

ksene=[]
list1=[]
list2=[]
for line in open("PHDISP.OUT"):
    i_line=line.split()
    if len(i_line):
       list1.append(float(i_line[0]))
       list2.append(abs(float(i_line[1])*2.194746313705e5))
    else:
       ksene.append([list1,list2])
       list1=[]
       list2=[]

#-------------------------------------------------------------------------------
# Read info about x-ticks position

bandlines=[]
fid=open("PHDLINES.OUT")
while 1:
    line=fid.readline()
    if not line:
        break
    i_line=line.split()
    bandlines.append(float(i_line[0]))
    fid.readline()
    fid.readline()
  
#-------------------------------------------------------------------------------
# Read info about x-ticks labels
  
ifile = open("input.xml","r")
iroot = et.parse(ifile).getroot()
ipath = iroot[3][0][0][0]
llist = []
for i in range(len(ipath)):
    label = " "
    if 'label' in ipath[i].attrib: label = ipath[i].attrib['label']
    if ( label == 'Gamma'): label = u'\u0393'
    llist.append(label)

#-------------------------------------------------------------------------------
# Settings for the plot

fontlabel=20
fontlegend=16
fonttick=16
fonttickx=20

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

ax1.xaxis.grid(True,which='major',color='k',linestyle='-',linewidth=2)
ax1.xaxis.set_label_position('bottom')
ax1.set_xticks(bandlines)

labels = ax1.set_xticklabels(llist)

ax1.text(-0.16,0.5,'Frequency [cm$\mathregular{^{-1}}$]',size=fontlabel,
        transform=ax1.transAxes,ha='center',va='center',rotation=90)

#-------------------------------------------------------------------------------
# Tick size

for line in ax1.get_xticklines() + ax1.get_yticklines(): line.set_markersize(9)
for line in ax1.get_yticklines(): 
    line.set_markeredgewidth(2)
    
ax1.yaxis.set_major_locator(MaxNLocator(7))

#-------------------------------------------------------------------------------
# Plot

ax1.plot(ksene[4][0],ksene[4][1],'r-',label="")
ax1.plot(ksene[3][0],ksene[3][1],'r-',label="")
ax1.plot(ksene[5][0],ksene[5][1],'r-',label="")

ax1.plot(ksene[0][0],ksene[0][1],'r-',label="")
ax1.plot(ksene[1][0],ksene[1][1],'r-',label="")
ax1.plot(ksene[2][0],ksene[2][1],'r-',label="")

ax1.yaxis.set_major_formatter(yfmt)

if ( narg > 0): plt.ylim(ymin=yymin)
if ( narg > 1): plt.ylim(ymax=yymax)

plt.xticks(size=fonttickx)
plt.yticks(size=fonttick)

ax1.text(1,1.05,plabel,size=fontlabel, color=mcolor,
             transform=ax1.transAxes,ha='right',va='center',rotation=0)

fig.subplots_adjust(left=None, bottom=None, right=None, wspace=None, hspace=None)

fig.savefig('PLOT.png',format='png',bbox_inches='tight',dpi=300)
fig.savefig('PLOT.pdf',format='pdf',bbox_inches=0)

sys.exit()    

#-------------------------------------------------------------------------------
