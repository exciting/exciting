#!/usr/bin/env python

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

#-------------------------------------------------------------------------------

def shell_value(variable,vlist,default):
    v = default
    e = False
    for i in range(len(vlist)):
        if ( vlist[i] == variable ): v = os.environ[variable] ; e = True ; break
    return v, e
    
#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])

#-------------------------------------------------------------------------------
  
###################
# Read input data #
###################

ha2ev = 27.2114

# KS
ksene=[]
list1=[]
list2=[]
for line in open("TDOS.OUT"):
    i_line=line.split()
    if len(i_line):
       list1.append(float(i_line[0])*ha2ev)
       list2.append(float(i_line[1])/ha2ev) # convert to eV
ksene.append([list1,list2])
del list1, list2

# G0W0
gwene=[]
list1=[]
list2=[]
for line in open("TDOS-QP.OUT"):
    i_line=line.split()
    if len(i_line):
       list1.append(float(i_line[0])*ha2ev)
       list2.append(float(i_line[1])/ha2ev) # convert to eV
gwene.append([list1,list2])
del list1, list2

################################################################################
################################################################################
################################################################################
    
figcolor = 'white'
dpi = 300
fig = plt.figure(figsize=(15,10),dpi=dpi)
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

ax1 = fig.add_axes([0.1,0.1,0.8,0.8])
# ax1.xaxis.grid(True,which='major',color='k',linestyle='-',linewidth=2)
ax1.set_xlabel('Energy (eV)')
ax1.set_ylabel(r'DOS (states/eV/bohr$^3$)')

for line in ax1.get_xticklines() + ax1.get_yticklines():
    line.set_markersize(10)
    line.set_markeredgewidth(2)

# Legend
ax1.plot(ksene[0][0],ksene[0][1],'b-',lw=3.0,label='KS')
ax1.plot(gwene[0][0],gwene[0][1],'r-',lw=3.0,label=r'G$_0$W$_0$')
leg=ax1.legend(bbox_to_anchor=(0.08,0.9),loc=2,borderaxespad=0.)
leg.draw_frame(False)

# add zero level
ymax = float(map(max,ksene[:][0])[1])
x0=[0.0,0.0]
y0=[0.0,ymax]
ax1.plot(x0,y0,'k:',lw=1.0)

# ax1.set_xlim(xmin,xmax)
# ax1.set_ylim(ymin,ymax)

# save the figure
plt.savefig('PLOT-dos.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT-dos.png', orientation='portrait',format='png',dpi=dpipng)

sys.exit()    
