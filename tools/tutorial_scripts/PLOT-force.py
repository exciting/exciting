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
import glob
import sys
import os

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

#-------------------------------------------------------------------------------

def sortstrain(s,e):
    ss=[]
    ee=[]
    ww=[]
    for i in range(len(s)): ww.append(s[i])
    ww.sort()
    for i in range(len(s)):
        ss.append(s[s.index(ww[i])])
        ee.append(e[s.index(ww[i])])
    return ss, ee

#-------------------------------------------------------------------------------

def shell_value(variable,vlist,default):
    v = default
    e = False
    for i in range(len(vlist)):
        if ( vlist[i] == variable ): v = os.environ[variable] ; e = True ; break
    return v, e
    
#-------------------------------------------------------------------------------

def leggi(filin,idf,a):
    x = [] ; f = True
    if (str(os.path.exists(filin))=='False'): 
        sys.exit("\nERROR: file "+filin+" not found!\n")
        return x, False
    os.system("grep -A"+a+" \""+idf+"\" "+str(filin)+" | grep \"at\" | grep \" "+a+" \" | tail -n1 > tempfile") 
    ifile = open("tempfile","r")
    x = ifile.readline().strip().split()[4:7]
    ifile.close()
    os.system("rm -f tempfile")
    return x, f
    
#-------------------------------------------------------------------------------

def leggishort(filin):
    f = open(filin,"r")
    x = float(f.readline().strip().split()[0])
    f.close()
    return x

#-------------------------------------------------------------------------------

ylabel = r'Force (cartesian) [Ha/Bohr]'
xlabel = u'Displacement $u$ [alat]'

#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

print "\n**Usage**:    PLOT-force.py [ATOM YMIN YMAX]\n"

#-------------------------------------------------------------------------------

a = str(1)
if (len(sys.argv) > 1): a = str(sys.argv[1])

#-------------------------------------------------------------------------------

label = "displ-"
if ( os.path.exists("strain-01") ): label = "strain-"
if ( os.path.exists("volume-01") ): label = "volume-"
if ( os.path.exists("alat-01") ): label = "alat-"

#-------------------------------------------------------------------------------

directoryroot = 'rundir-'
list_dir = sorted(glob.glob(directoryroot+"*"))

xx = [] ; y1 = [] ; y2 = [] ; y3 = [] 

for idir in range(len(list_dir)):    
    c = [] 
    fileinp = list_dir[idir]+'/INFO.OUT'
    c, fcheck = leggi(fileinp,"Total atomic forces",a)
    if (fcheck): 
        y1.append(float(c[0]))
        y2.append(float(c[1]))
        y3.append(float(c[2]))
        xx.append(leggishort(label+list_dir[idir][-2:]))

if ( label == 'displ-' ):     
    of = open("force-vs-displacement","w")
    fmt = '%16.8f'
    for i in range(len(xx)):
        print>>of, (fmt%xx[i]), (fmt%y1[i]), (fmt%y2[i]), (fmt%y3[i])
    of.close()
        
#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])

#-------------------------------------------------------------------------------
# manipulate data for a better plot

xmin = min(xx)
xmax = max(xx)

ymin = min(y1+y2+y3)
ymax = max(y1+y2+y3)

dxx  = abs(xmax-xmin)/18
dyy  = abs(ymax-ymin)/15

xmin = xmin-dxx
xmax = xmax+dxx

ymin = ymin-dyy
ymax = ymax+dyy

#-------------------------------------------------------------------------------
# set defauls parameters for the plot

fontlabel=20
fonttick=16
fonttext=14

params = {'ytick.minor.size': 6,
          'xtick.major.pad': 8,
          'ytick.major.pad': 4,
          'patch.linewidth': 2.,
          'axes.linewidth': 2.,
          'lines.linewidth': 1.8,
          'lines.markersize': 8.0,
          'axes.formatter.limits': (-4, 6)}

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
 
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(6)
    line.set_markeredgewidth(2)      

plt.xticks(size=fonttick)
plt.yticks(size=fonttick)
pyl.grid(True)

plt.plot(xx,y1,'ro-',label=u'$F_x$')
plt.plot(xx,y2,'bo-',label=u'$F_y$')
plt.plot(xx,y3,'go-',label=u'$F_z$')

plt.plot(xx,y3,'go-')
plt.plot(xx,y2,'bo-')
plt.plot(xx,y1,'ro-')

#-------------------------------------------------------------------------------

plt.legend(bbox_to_anchor=(1.03, 1), loc=2, borderaxespad=0., numpoints=1)

ax.yaxis.set_major_formatter(yfmt)

if (abs(ymax-ymin) < 0.000000001): 
    ymax=ymax+0.1
    ymin=ymin-0.1

if (len(sys.argv) > 2): ymin = float(sys.argv[2])
if (len(sys.argv) > 3): ymax = float(sys.argv[3])

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

ax.xaxis.set_major_locator(MaxNLocator(7))

ax.set_axisbelow(True) 

plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------







