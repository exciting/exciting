#!/usr/bin/env python

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt

###################
# Read input data #
###################


# KS states from BAND.OUT
ksene=[]
list1=[]
list2=[]
for line in open("BAND.OUT"):
    i_line=line.split()
    if len(i_line):
       list1.append(float(i_line[0]))
       list2.append(float(i_line[1])) # convert to eV
    else:
       ksene.append([list1,list2])
       list1=[]
       list2=[]

# G0W0 states from BAND-QP.OUT
gwene=[]
list1=[]
list2=[]
for line in open("BAND-QP.OUT"):
    i_line=line.split()
    if len(i_line):
       list1.append(float(i_line[0]))
       list2.append(float(i_line[1])) # convert to eV
    else:
       gwene.append([list1,list2])
       list1=[]
       list2=[]


# Read info about x-ticks position
bandlines=[]
fid=open("BANDLINES.OUT")
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

# Total DOS for KS and GW states from TDOS-KS-QP.OUT
ene=[]
ksdos=[]
gwdos=[]
for line in open("TDOS-KS-QP.OUT"):
    i_line=line.split()
    ene.append(float(i_line[0]))
    ksdos.append(float(i_line[1]))
    gwdos.append(float(i_line[2]))


################################################################################
################################################################################
################################################################################
    
figcolor = 'white'
dpi = 100
fig = plt.figure(figsize=(20,10),dpi=dpi)
fig.figurePatch.set_edgecolor(figcolor)
fig.figurePatch.set_facecolor(figcolor)

mpl.rcParams['axes.linewidth'] = 3.0 # set the value globally
mpl.rcParams['grid.linewidth'] = 1.5
mpl.rcParams['xtick.labelsize'] = 25
mpl.rcParams['ytick.labelsize'] = 25
mpl.rcParams['axes.edgecolor'] = 'black'
mpl.rcParams['axes.labelsize'] = '30'     # fontsize of the x any y labels
mpl.rcParams['axes.labelcolor'] = 'black'
mpl.rcParams['axes.axisbelow'] = 'True'   # whether axis gridlines and ticks are below
                                          # the axes elements (lines, text, etc)
mpl.rcParams['legend.fontsize'] = '25'
mpl.rcParams['legend.fontsize'] = '25'

plt.rcParams['xtick.major.pad'] = '10'

#fig.text(0.84,0.67,r'$[1\bar{1}0]$',fontsize=35)

#############################
##    Bandstructure plot   ##
#############################

ax1 = fig.add_axes([0.1,0.1,0.6,0.8])
ax1.xaxis.grid(True,which='major',color='k',linestyle='-',linewidth=2)
ax1.set_xlabel('Bandstructure')
ax1.xaxis.set_label_position('top')
ax1.set_xticks(bandlines)
labels = ax1.set_xticklabels(('W','L',r'$\Gamma$','X','W','K'))
ax1.set_ylabel('Energy [Ha]')

#for line in ax1.get_xticklines() + ax1.get_yticklines():
#    line.set_markersize(10)
#    line.set_markeredgewidth(2)

ymin=1000.0
ymax=-1000.0
for i in range(len(gwene)):
    y=min(gwene[i][1])
    if (y<ymin):
        ymin=y
    y=max(gwene[i][1])
    if (y>ymax):
        ymax=y

for i in range(len(gwene)):
    ax1.plot(ksene[i][0],ksene[i][1],'b.-',lw=2.0)
    ax1.plot(gwene[i][0],gwene[i][1],'r.-',lw=2.0)

# add zero level
x0=[0.0,max(ksene[0][0])]
y0=[0.0,0.0]
ax1.plot(x0,y0,'k:',lw=1.0)

ax1.set_xlim(0,max(ksene[0][0]))
ax1.set_ylim(ymin,ymax)

#############################
##         DOS plot        ##
#############################

rect = [0.72, 0.1, 0.20, 0.8]
ax2 = fig.add_axes(rect)
labels = ax2.set_xticklabels((''))
ax2.set_xlabel('DOS')
ax2.xaxis.set_label_position('top')
ax2.set_ylabel('Energy [Ha]')
# Tick size
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position('right')

for line in ax2.get_xticklines() + ax2.get_yticklines():
    line.set_markersize(10)
    line.set_markeredgewidth(2)

# Legend

ax2.plot(ksdos,ene,'b.-',lw=2.0,label='KS')
ax2.plot(gwdos,ene,'r.-',lw=2.0,label='GW')

leg=ax2.legend(bbox_to_anchor=(0.5,0.2),loc=2,borderaxespad=0.)
leg.draw_frame(False)

# add zero level
x0=[0.0,0.5*max(ksdos)]
y0=[0.0,0.0]
ax2.plot(x0,y0,'k:',lw=1.0)

xmin=0.0
xmax=0.5*max(ksdos)
ax2.set_xlim(xmin,xmax)
ax2.set_ylim(ymin,ymax)


fig.savefig('Figure.png',format='png',bbox_inches=0,dpi=50)

plt.show()
sys.exit()    
