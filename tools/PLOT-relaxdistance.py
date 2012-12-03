#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   sys   import stdin
from   math  import sqrt
from   math  import factorial
from   pylab import *
import matplotlib.transforms as ptf
import matplotlib.ticker     as ptk 
import matplotlib.pyplot     as plt
import pylab                 as pyl
import numpy
import sys
import os

#-------------------------------------------------------------------------------

def shell_value(variable,vlist,default):
    v = default
    e = False
    for i in range(len(vlist)):
        if ( vlist[i] == variable ): v = os.environ[variable] ; e = True ; break
    return v, e

#-------------------------------------------------------------------------------

xlabel  = u'Optimization steps'
ylabel  = r'Relative coordinate [crystal]'

#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])

#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

if (narg < 1): 
    print "\nIncorrect number of arguments. **Usage**:\n\n",
    print "PLOT-relaxdistance.py DIRECTORYNAME [YMIN YMAX]\n"
    sys.exit()

label = str(sys.argv[1])

#-------------------------------------------------------------------------------

inpf    = current+"/"+rlabel+label+'/INFO.OUT'
if (label == 'r'): inpf=rundir+'/xc-rundir/INFO.OUT'
     
if (str(os.path.exists(inpf))=='False'): 
    sys.exit("\nERROR: file "+inpf+" not found!\n")
 
#-------------------------------------------------------------------------------
   
ylimits = []
for i in range(2,len(sys.argv)): ylimits.append(float(sys.argv[i]))
    
#-------------------------------------------------------------------------------

atom1=1
atom2=2

os.system("grep \"   "+str(atom1)+" : \" "+str(inpf)+" > tempfile1")
ifile1 = open("tempfile1","r")
os.system("grep \"   "+str(atom2)+" : \" "+str(inpf)+" > tempfile2")
ifile2 = open("tempfile2","r")

#-------------------------------------------------------------------------------
# set defauls parameters for the plot

fontlabel=20
fonttick=16

params = {'ytick.minor.size': 6,
          'xtick.major.pad': 8,
          'ytick.major.pad': 4,
          'patch.linewidth': 2.,
          'axes.linewidth': 2.,
          'lines.linewidth': 1.8,
          'lines.markersize': 8.0,
          'axes.formatter.limits': (-4, 4)}
 
plt.rcParams.update(params)

plt.subplots_adjust(left=0.22, right=0.78,
                    bottom=0.18, top=0.88,
                    wspace=None, hspace=None)
                    
yfmt = ptk.ScalarFormatter(useOffset=True,useMathText=True)

fig  = matplotlib.pyplot.figure(1, figsize=(8,5.5)) 

ax   = fig.add_subplot(111)

ax.text(0.5,-0.13,xlabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=0)
ax.text(-0.23,0.5,ylabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=90)

#-------------------------------------------------------------------------------

x = [] ; d1 = [] ; d2 = [] ; d3 = [] 

iter=0

while True:
    line1 = ifile1.readline().strip()
    line2 = ifile2.readline().strip()   
    if len(line1) == 0: break
    iter+=1
    d1.append(float(line2.split()[2])-float(line1.split()[2]))
    d2.append(float(line2.split()[3])-float(line1.split()[3]))
    d3.append(float(line2.split()[4])-float(line1.split()[4]))
    x.append(float(iter))

xmin = 1-iter/20. ; xmax = iter+iter/20.

#-------------------------------------------------------------------------------
# manipulate data for a better plot

amax  = max(max(d1),max(d2),max(d3))
amin  = max(min(d1),min(d2),min(d3))
aaaa  = max(abs(amax),abs(amin))

newl  = True
rmin  = d1[0] 
if (aaaa < 1000): newl = False ; rmin = 0. 

srmin = u'\u2013 '+str(abs(rmin))
if (rmin > 0): srmin = u'+ '+str(rmin)

for i in range(len(d1)): d1[i]=(d1[i]-rmin)
for i in range(len(d2)): d2[i]=(d2[i]-rmin)
for i in range(len(d3)): d3[i]=(d3[i]-rmin)

ymin = min(min(d1),min(d2),min(d3))
ymax = max(max(d1),max(d2),max(d3))

if (len(ylimits) == 1): 
    ymin = float(ylimits[0])-rmin
if (len(ylimits) > 1): 
    ymin = float(ylimits[0])-rmin ; ymax = float(ylimits[1])-rmin 
    
dyy  = abs(ymax-ymin)/18
ymin = ymin-dyy ; ymax = ymax+dyy

#-------------------------------------------------------------------------------

if (newl): ax.text(0.11,1.03,srmin,size=fonttext,
        transform=ax.transAxes,ha='left',va='center',rotation=0)
 
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(6)
    line.set_markeredgewidth(2)
 
plt.xticks(size=fonttick)
plt.yticks(size=fonttick)
pyl.grid(True)

#-------------------------------------------------------------------------------

os.system("rm tempfile1") ; os.system("rm tempfile2")

plt.plot(x,d3,'go-',label=u'$\Delta$3')
plt.plot(x,d2,'bo-',label=u'$\Delta$2')
plt.plot(x,d1,'ro-',label=u'$\Delta$1')

#-------------------------------------------------------------------------------

iloc = 1
mdum = min(d1[-1],d2[-1],d3[-1])
ydum = (ymax-ymin)/2.2+ymin
if (mdum > ydum): iloc = 4
plt.legend(loc=iloc,borderaxespad=.8)

plt.legend(bbox_to_anchor=(1.03, 1), loc=2, borderaxespad=0.)

ax.yaxis.set_major_formatter(yfmt)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

ax.set_axisbelow(True) 

plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------





