#!/usr/bin/python
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from lxml import etree
import os
import sys

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
fname="MOKE_NLF_FXCRPA_QMT001.OUT"
#-------------------------------------------------------------------------------
#Parse LOSS function data files
xdata=[]
ydata=[]
labels=[]
legends=[]
function=1

print "Parsing "+fname

tree=etree.parse(fname)

for elem in tree.xpath('/%s/map'%(rootelement)):
    xdata[i].append(float(elem.attrib["variable1"]))
    ydata[i].append(float(elem.attrib["function%d"%(function)]))

#-------------------------------------------------------------------------------
#Plot LOSS function/s 
colors=['k','r','g','b','y','c','m']
fig=plt.figure(1,figsize=(8,5.5))

params = {'font.size':15,
          'xtick.major.size': 5,
          'ytick.major.size': 5,
          'patch.linewidth': 1.5,
          'axes.linewidth': 2.,
          'axes.formatter.limits': (-4, 6),
          'lines.linewidth': 1.8,
          'lines.markeredgewidth':2.0,
          'lines.markersize':18,
          'legend.fontsize':11,
          'legend.borderaxespad':1,
          'legend.borderpad':0.5,
          'savefig.dpi':80}

plt.rcParams.update(params)

ax=fig.add_subplot(111)

for i in range(nfiles):
    ax.plot(xdata[i],ydata[i],colors[np.mod(i,7)],label=legends[i])

ax.legend(loc=2)
#ax.legend()

ax.set_xlim(0.0,54.0)
ax.set_ylim(0)
ax.set_xlabel(str.capitalize(labels[0]["xlabel"])+" [eV]")
ax.set_ylabel(str.capitalize(labels[0]["ylabel"]))

plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
