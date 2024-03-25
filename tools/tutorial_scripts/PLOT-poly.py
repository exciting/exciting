#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   sys   import stdin
from   math  import sqrt
from   math  import factorial
from   pylab import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ptk 
import pylab             as pyl
import numpy
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

if (str(os.path.exists('energy-vs-volume'))=='False'): 
    sys.exit("ERROR: file energy-vs-volume not found!\n")

#-------------------------------------------------------------------------------

order_of_fit = 6

narg  = len(sys.argv)-1

if (narg>0): order_of_fit = int(sys.argv[1])

#-------------------------------------------------------------------------------

print

#print "==============================="
#print "Lattice symmetry codes"
#print "-------------------------------"
#print "1 --> Simple cubic (sc)"
#print "2 --> Body-centered cubic (bcc)"
#print "3 --> Face-centered cubic (fcc)"
#print "-------------------------------"
#print "0 --> Others"
#print "===============================\n"

scheck = "3"\
     #raw_input("Enter lattice symmetry code [default 0] >>>> ").replace(" ", "") 

isym   = 0
factor = 1
if ( scheck == "1" ): isym = 1 ; factor=1 ; slabel = "(sc) "
if ( scheck == "2" ): isym = 2 ; factor=2 ; slabel = "(bcc)"
if ( scheck == "3" ): isym = 3 ; factor=4 ; slabel = "(fcc)"
#print "Verification lattice symmetry code      >>>>", isym

#-------------------------------------------------------------------------------

energy = []
strain = []

#-------------------------------------------------------------------------------

input_energy = open('energy-vs-volume',"r")

while True:
    line = input_energy.readline()
    line = line.strip()
    if len(line) == 0: break
    energy.append(float(line.split()[1]))
    strain.append(float(line.split()[0]))

strain,energy=sortstrain(strain,energy)

#-------------------------------------------------------------------------------

bohr_radius     = 0.529177
joule2hartree   = 4.3597482
joule2rydberg   = joule2hartree/2.
unitconv        = joule2hartree/bohr_radius**3*10.**3

#-------------------------------------------------------------------------------

fitr = numpy.polyfit(strain,energy,order_of_fit)
curv = numpy.poly1d(fitr)
bulk = numpy.poly1d(numpy.polyder(fitr,2))
bpri = numpy.poly1d(numpy.polyder(fitr,3))
vmin = numpy.roots(numpy.polyder(fitr))

dmin=[]
for i in range(len(vmin)):
    if (abs(vmin[i].imag) < 1.e-10): 
        if (strain[0] <= vmin[i] and vmin[i] <= strain[-1]): 
            if(bulk(vmin[i]) > 0): dmin.append(vmin[i].real)

xvol = numpy.linspace(strain[0],strain[-1],100)

chi = 0
for i in range(len(energy)): 
    chi=chi+(energy[i]-curv(strain[i]))**2
chi=sqrt(chi)/len(energy)

#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])

#-------------------------------------------------------------------------------

xlabel = u'Volume [Bohr\u00B3]'
ylabel = r'Energy [Ha]'
if (os.path.exists('quantum-espresso')): ylabel = r'Energy [Ry]'
if (os.path.exists('vasp')): ylabel = r'Energy [Ry]'

#-------------------------------------------------------------------------------

fontlabel=20
fonttick=16

params = {'ytick.minor.size': 6,
          'xtick.major.pad': 8,
          'ytick.major.pad': 4,
          'patch.linewidth': 2.,
          'axes.linewidth': 2.,
          'lines.linewidth': 1.8,
          'lines.markersize': 8.0,
          'axes.formatter.limits': (-4, 6)}

plt.rcParams.update(params)
plt.subplots_adjust(left=0.21, right=0.93,
                    bottom=0.18, top=0.88,
                    wspace=None, hspace=None)
                           
yfmt = ptk.ScalarFormatter(useOffset=True,useMathText=True)                           
                           
figure = plt.figure(1, figsize=(8,5.5))  
ax     = figure.add_subplot(111)
ax.text(0.5,-0.13,xlabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center')
ax.text(-0.19,0.5,ylabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=90)
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(6)
    line.set_markeredgewidth(2)
plt.xticks(size=fonttick)
plt.yticks(size=fonttick)
pyl.grid(True)
plt.plot(xvol,curv(xvol),'b-',label='n='+str(order_of_fit)+' fit')
plt.plot(strain,energy,'go',label='calculated')
plt.plot(dmin,curv(dmin),'ro')
plt.legend(loc=9,borderaxespad=.8,numpoints=1)

ymax  = max(max(curv(xvol)),max(energy))
ymin  = min(min(curv(xvol)),min(energy))
dxx   = abs(max(xvol)-min(xvol))/18
dyy   = abs(ymax-ymin)/18
ax.yaxis.set_major_formatter(yfmt)
ax.set_xlim(min(xvol)-dxx,max(xvol)+dxx)
ax.set_ylim(ymin-dyy,ymax+dyy)

ax.xaxis.set_major_locator(MaxNLocator(7))

plt.savefig('PLOT.ps', orientation='portrait',format='eps')
plt.savefig('PLOT.png',orientation='portrait',format='png',dpi=dpipng)

#-------------------------------------------------------------------------------

#print 
#print "##############################################\n"
if (len(dmin) > 1): 
    print 
    print "##############################################\n"
    print "WARNING: Multiple minima are found!\n"
    print "##############################################\n"

fmt='%11.5f'
amt='%10.4f'
bmt='%9.3f'
pmt='%16.10f'
lmt='%10.2f'

for i in range(len(dmin)):
    v0=dmin[len(dmin)-1-i]
    a0sc=(1*v0)**(0.33333333333)
    abcc=(2*v0)**(0.33333333333)
    afcc=(4*v0)**(0.33333333333)
    a0=(factor*v0)**(0.33333333333)
    b0=bulk(v0)*v0*unitconv
    bp=-(1+v0*bpri(v0)/bulk(v0))
#    print 'Optimal volume   = ', fmt%(v0), '[Bohr^3]'
#    if (isym > 0): print 'Lattice constant =', slabel, afmt%(a0), '[Bohr]'
#    print 'Bulk modulus     = ', fmt%(b0), '[GPa]'
#    print
#    print 'Log(chi)         = ', lmt%(log10(chi))
#    print
    print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print "     V0        B0         Bp        a-sc       a-bcc      a-fcc     log(chi)"
    print fmt%(v0), bmt%(b0), bmt%(bp),  
    print amt%(a0sc), amt%(abcc), amt%(afcc), lmt%(log10(chi))
    print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print
    
if ( len(dmin) == 0): 
    print
    print "WARNING: No minimum in the given xrange!\n"
    print "##############################################\n"
 
if (showpyplot): plt.show()
#-------------------------------------------------------------------------------












