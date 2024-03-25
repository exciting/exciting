#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   pylab import *
import sys
import matplotlib        as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ptk 
import pylab             as pyl
import numpy
import os

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

#-------------------------------------------------------------------------------
    
def change(n1,n2,i1,i2,x):
    step = 1
    if (x != 0): step = -1
    for i in range(i1,i2+1,step):
        k1 = i
        k2 = i
        if (x != 0): 
            k1 = len(ksene[n1][1])-i
            k2 = len(ksene[n2][1])-i     
       #if (abs(ksene[n1][0][k1]-ksene[n2][0][k2]) > 0.0001):
       #    print n1,k1,ksene[n1][0][k1],"  ",n2,k2, ksene[n2][0][k2]
        adum = ksene[n1][1][k1]
        ksene[n1][1][k1] = ksene[n2][1][k2]
        ksene[n2][1][k2] = adum
        adum = ksene[n1][0][k1]
        ksene[n1][0][k1] = ksene[n2][0][k2]
        ksene[n2][0][k2] = adum
    return

def eraser(x,y,i0,i1):
    z = [] ; f = []
    z.append(x[i0-1])
    z.append(x[i1+1])
    f.append(y[i0-1])
    f.append(y[i1+1])
    b, a = numpy.polyfit(z,f,1)
    for i in range(i0,i1+1):
        y[i] = a+b*x[i]
    return
    
#-------------------------------------------------------------------------------
# Read data

ksene=[]
list1=[]
list2=[]
for line in open("PHDISP.OUT"):
    i_line=line.split()
    if len(i_line):
       list1.append(float(i_line[0]))
       list2.append(float(i_line[1]))
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
          'lines.linewidth': 2.5,
          'lines.markersize': 5.0,
          'axes.formatter.limits': (-5, 6)}

plt.rcParams.update(params)

plt.subplots_adjust(left=0.20, right=0.93,
                    bottom=0.18, top=0.92,
                    wspace=None, hspace=None)
                    
yfmt = ptk.ScalarFormatter(useOffset=True,useMathText=True)

fig  = matplotlib.pyplot.figure(1, figsize=(8,5.5)) 

ax1  = fig.add_subplot(111)

ax1.xaxis.grid(True,which='major',color='k',linestyle='-',linewidth=2)
ax1.xaxis.set_label_position('bottom')
ax1.set_xticks(bandlines)

labels = ax1.set_xticklabels((u'\u0393','K','X',u'\u0393','L'))
ax1.text(-0.16,0.5,r'Grueneisen parameter',size=fontlabel,
        transform=ax1.transAxes,ha='center',va='center',rotation=90)

#-------------------------------------------------------------------------------
# Tick size

for line in ax1.get_xticklines() + ax1.get_yticklines(): line.set_markersize(9)
for line in ax1.get_yticklines(): 
    line.set_markeredgewidth(2)
    
ax1.yaxis.set_major_locator(MaxNLocator(7))

#-------------------------------------------------------------------------------
# Eliminate crosses

change(2, 4, 258,  -1, 1)
eraser(ksene[2][0],ksene[2][1],237,243)
eraser(ksene[4][0],ksene[4][1],237,243)

change(1, 2, 141, 154, 0)
eraser(ksene[1][0],ksene[1][1],141,141)
eraser(ksene[2][0],ksene[2][1],141,141)
eraser(ksene[1][0],ksene[1][1],154,155)
eraser(ksene[2][0],ksene[2][1],154,155)

change(3, 4, 125, 167, 0)
eraser(ksene[3][0],ksene[3][1],124,133)
eraser(ksene[4][0],ksene[4][1],124,133)
eraser(ksene[3][0],ksene[3][1],165,171)
eraser(ksene[4][0],ksene[4][1],165,171)

change(4, 5,   0,  98, 0)
eraser(ksene[4][0],ksene[4][1],95,100)
eraser(ksene[5][0],ksene[5][1],95,100)

change(4, 2,   0,  214, 0)

#-------------------------------------------------------------------------------
# Plot

ax1.plot(ksene[2][0],ksene[2][1],'b-',label="")
ax1.plot(ksene[3][0],ksene[3][1],'b-',label="")
ax1.plot(ksene[5][0],ksene[5][1],'b-',label="")

ax1.plot(ksene[0][0],ksene[0][1],'r-',label="")
ax1.plot(ksene[1][0],ksene[1][1],'r-',label="")
ax1.plot(ksene[4][0],ksene[4][1],'r-',label="")

ax1.yaxis.set_major_formatter(yfmt)

ax1.set_xlim(0,max(ksene[0][0])+1.e-7)
#ax1.set_xlim(1.4,3.0688109351)
#ax1.set_xlim(1.32,1.36)
ax1.set_ylim(0,1.6)

plt.xticks(size=fonttickx)
plt.yticks(size=fonttick)

fig.subplots_adjust(left=None, bottom=None, right=None, wspace=None, hspace=None)

fig.savefig('PLOT.png',format='png',bbox_inches=0,dpi=300)
fig.savefig('PLOT.pdf',format='pdf',bbox_inches=0)

sys.exit()    

#-------------------------------------------------------------------------------
