#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   lxml  import etree
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

def leggi(filin,idf,a1,a2):
    os.system("grep -A"+a1+" \""+idf+"\" "+str(filin)+" | grep \"at\" | grep \" "+a1+" \" | tail -n1 > tempfile1") 
    ifile1 = open("tempfile1","r")
    os.system("grep -A"+a2+" \""+idf+"\" "+str(filin)+" | grep \"at\" | grep \" "+a2+" \" | tail -n1 > tempfile2") 
    ifile2 = open("tempfile2","r")
    x = [] ; y = []
    x = ifile1.readline().strip().split()[4:7]
    y = ifile2.readline().strip().split()[4:7]
    ifile1.close()
    ifile2.close()
    os.system("rm -f tempfile1 tempfile2")
    return x,y

#-------------------------------------------------------------------------------

def leggishort(filin):
    f = open(filin,"r")
    x = float(f.readline().strip().split()[0])
    f.close()
    return x

#-------------------------------------------------------------------------------

def relax(filin,idf):
    check = True
    os.system("grep \""+idf+"\" "+filin+" > tempfile")
    f = open("tempfile","r")
    x = f.readline().strip().split()
    if (len(x) < 1): check = False
    os.system("rm -f tempfile")
    return check
    
#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])

#-------------------------------------------------------------------------------

lespresso = os.path.exists("quantum-espresso")
if (lespresso): sys.exit("\n ERROR: Quantum ESPRESSO version not yet implemented!\n")

#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

print "\n**Usage**:    PLOT-optimized-geometry.py [ATOM1 ATOM2 YMIN YMAX]\n"

#-------------------------------------------------------------------------------

a1 = str(1)
if (len(sys.argv) > 1): a1 = str(sys.argv[1])
a2 = str(2)
if (len(sys.argv) > 2): a2 = str(sys.argv[2])

#-------------------------------------------------------------------------------

list_dir = sorted(glob.glob('rundir-*'))

#-------------------------------------------------------------------------------

if (str(os.path.exists(current+'/'+list_dir[0]+'/input.xml'))=='False'): 
    sys.exit("ERROR: Input file "+current+"/"+list_dir[0]+"/input.xml not found!\n")

acoord = "lattice"

input_obj = open(current+"/"+list_dir[0]+"/input.xml","r")
input_doc = etree.parse(input_obj)
input_rut = input_doc.getroot()
 
xml_cartesian = map(str,input_doc.xpath('/input/structure/@cartesian'))
if (xml_cartesian == []):
    acoord = "lattice"
else:
    if (xml_cartesian[0] == "true"): 
        acoord = "cartesian"   

#-------------------------------------------------------------------------------

xx = [] ; x0 = []
y1 = [] ; y2 = [] ; y3 = [] 
z1 = [] ; z2 = [] ; z3 = [] 

for idir in range(len(list_dir)):
    r1 = [] ; r2 = []    
    c1 = [] ; c2 = []
    fileinp = list_dir[idir]+'/INFO.OUT'
    r1,r2 = leggi(fileinp,"Atomic positions (",a1,a2)
    check = relax(fileinp,"Atomic positions at this step")
    if (check):
        c1,c2 = leggi(fileinp,"Atomic positions at this step",a1,a2)
    else:
        c1,c2 = leggi(fileinp,"Atomic positions (",a1,a2)
    
    y1.append(float(c2[0])-float(c1[0]))
    y2.append(float(c2[1])-float(c1[1]))
    y3.append(float(c2[2])-float(c1[2]))
    xx.append(leggishort('strain-'+list_dir[idir][-2:]))
    if (idir == 0 or idir == len(list_dir)-1):
        z1.append(float(r2[0])-float(r1[0]))
        z2.append(float(r2[1])-float(r1[1]))
        z3.append(float(r2[2])-float(r1[2]))
        x0.append(leggishort('strain-'+list_dir[idir][-2:]))
    
ylabel  = r'Relative coordinate ('+acoord+')'
xlabel  = r'Lagrangian strain'

#-------------------------------------------------------------------------------
# manipulate data for a better plot

xmin = min(xx+x0)
xmax = max(xx+x0)

ymin = min(y1+y2+y3+z1+z2+z3)
ymax = max(y1+y2+y3+z1+z2+z3)

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

plt.plot(xx,y1,'ro--',label=u'$\Delta$1')
plt.plot(xx,y2,'bs--',label=u'$\Delta$2')
plt.plot(xx,y3,'gd--',label=u'$\Delta$3')

plt.plot(x0,z1,'r-',label=u'$\Delta$1$_{ref}$')
plt.plot(x0,z2,'b-',label=u'$\Delta$2$_{ref}$')
plt.plot(x0,z3,'g-',label=u'$\Delta$3$_{ref}$')

plt.plot(x0,z3,'g-')
plt.plot(x0,z2,'b-')
plt.plot(x0,z1,'r-')

plt.plot(xx,y3,'gd--')
plt.plot(xx,y2,'bs--')
plt.plot(xx,y1,'ro--')

geo_output = open("opt-relative-geometry-lattice","w")
for i in range(len(xx)): 
    print >>geo_output, xx[i], 
    print >>geo_output, y1[i]-z1[0], 
    print >>geo_output, y2[i]-z2[0], 
    print >>geo_output, y3[i]-z3[0] 
geo_output.close()

#-------------------------------------------------------------------------------

#plt.legend(loc=9,borderaxespad=.8,numpoints=1)
plt.legend(bbox_to_anchor=(1.03, 1), loc=2, borderaxespad=0., numpoints=1)

ax.yaxis.set_major_formatter(yfmt)

#ylimits = []
#for i in range(1,len(sys.argv)): ylimits.append(float(sys.argv[i]))

#if (len(ylimits) == 1): ymin = float(ylimits[0])
#if (len(ylimits) > 1): ymin = float(ylimits[0]); ymax = float(ylimits[1]) 

if (abs(ymax-ymin) < 0.000000001): 
    ymax=ymax+0.1
    ymin=ymin-0.1

if (len(sys.argv) > 3): ymin = float(sys.argv[3])
if (len(sys.argv) > 4): ymax = float(sys.argv[4])

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

ax.xaxis.set_major_locator(MaxNLocator(7))

ax.set_axisbelow(True) 

plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------





