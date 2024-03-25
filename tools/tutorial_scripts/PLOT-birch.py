#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   lxml  import etree
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
import math 
import numpy
from   scipy import *
from   scipy.optimize import leastsq

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

#-------------------------------------------------------------------------------

def inifit(x,y):
    p=[] ; a,b,c = polyfit(x,y,2)
    p.append(-b/(2*a))               # v0
    p.append(a*p[0]**2 + b*p[0] + c) # e0
    p.append((2*a*p[0]))             # b0
    p.append(2.)                     # bp
    return p

def bmeos(v,p):
    v0 = p[0] ; e0 = p[1] ; b0 = p[2] ; bp = p[3]  
    vv = (v0/v)**(2./3.)
    ff = vv - 1.
    ee = e0 + 9.*b0*v0/16. * ( ff**3*bp + ( 6.-4.*vv )*ff**2 )
    return ee    

def residuals(p,e,v):
    return e - bmeos(v,p)

def bmfit(x,y):
    p = [0,0,0,0] ; p = inifit(x,y)
    return leastsq(residuals, p, args=(y, x))

#-------------------------------------------------------------------------------

def printderiv(damax,bb,nfit,output_file):
    print >> output_file, '%16.10f'%(damax),
    for j in range(0,nfit):
        if (bb[j] != 99999999):
            print >> output_file, '%14.6f'%(bb[j]),
    print >> output_file
    return 

def printnice(etam,bb,nfit,nderiv,orderstep,strainpoints,dc,dl):
    print
    print "###########################################\n"
    print "Fit data-----------------------------------\n"
    print "Deformation code             ==>", dc#'%2i'%(dc)
    print "Deformation label            ==>", dl
    print "Maximum value of the strain  ==>", '%10.8f'%(etam) 
    print "Number of strain values used ==>", '%2i'%(strainpoints),"\n" 
    print "Fit results for the derivative of order", '%3i'%(nderiv),"\n"  
    for j in range(0,nfit):
        order=nderiv+j*orderstep
        if (bb[j] != 99999999):
            print "Polynomial of order", '%2i'%(order), "==>", 
            print '%10.2f'%(bb[j]), "[GPa]"
    print
    return 

def sortstrain(s,e):
    ss=[] ; ee=[] ; ww=[]
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
    
def eos(v,v0,e0,b0,bp):
    vv = (v0/v)**(2./3)
    ff = vv - 1.
    ee = e0 + 9.*b0*v0/16. * ( ff**3*bp + ( 6.-4.*vv )*ff**2 )
    return ee   
    
def pickfromfile(filename,nrow,ncol):
    pfile = open(filename,"r")
    for i in range(nrow): line = pfile.readline().strip()
    x = line.split()[ncol-1]
    pfile.close()
    return x 
    
def vol2eta(v0,vzero,n):
    v = v0
    if (n == 0): return v
    v = v0/vzero
    y = 1./float(n)
    e = v**y -1.0
    x = e + e**2/2.0
    return x
    
#-------------------------------------------------------------------------------

maxstrain = 1.
args = sys.argv[1:]
if ( (len(args) > 0) and ("-maxstrain" in args) ): 
    ind = args.index("-maxstrain")
    maxstrain = float(args[ind+1])
    
#-------------------------------------------------------------------------------
   
lexciting = os.path.exists('exciting')
lespresso = os.path.exists('quantum-espresso')
lvasp     = os.path.exists('vasp')
lrydberg  = lespresso or lvasp 
lplanar   = os.path.exists('planar')
lexsinfo  = os.path.exists('INFO-elastic-constants')
lexsev    = os.path.exists('energy-vs-volume')
lexses    = os.path.exists('energy-vs-strain')
lsource   = os.path.exists('source.xml')

if (not(lexsev) and not(lexses)): 
    sys.exit("\n ERROR: file energy-vs-volume or energy-vs-strain not found!\n")
    
#-------------------------------------------------------------------------------

factor    = 2.
if (lrydberg): factor = 1.

if (lexsev): 
    print "\n Input file is \"energy-vs-volume\"."
    input_energy = open('energy-vs-volume',"r")
    
else:
    print "\n Input file is \"energy-vs-strain\"."
    if (not(lexsinfo)): sys.exit("\n ERROR: file INFO-elastic-constants not found!\n") 
    input_energy = open('energy-vs-strain',"r")
    #________________________________________
    # check dimentionality of the deformation
    #
    deformation = int(pickfromfile("INFO-elastic-constants",5,4))
    if   (deformation == 0): dimension = 3
    elif (deformation == 8): dimension = 2
    elif ( (deformation == 1) or (deformation == 2) or (deformation == 3) ): dimension = 1
    else: sys.exit("\n ERROR: deformation type "+str(deformation)+" not allowed!\n")
    bfactor = dimension**2
    #___________________________
    # read volume at zero strain
    #
    vzero = float(pickfromfile("INFO-elastic-constants",4,7))
    
#-------------------------------------------------------------------------------

if (lexses and lplanar):    
    print " Modified version for strained planar systems!"
    #_____________________
    # read alat and covera
    #
    alat   = float(pickfromfile('planar',1,1))
    covera = float(pickfromfile('planar',1,2))
    factor = factor*alat*covera*5.2917721092e-2
    
#-------------------------------------------------------------------------------

bohr_radius     = 0.529177
joule2hartree   = 4.3597482
joule2rydberg   = joule2hartree/2.
unitconv        = joule2hartree/bohr_radius**3*10.**3*factor

electron_charge = 1.602176565e-19 
bohr_radius     = 5.2917721092e-11
rydberg2ev      = 13.605698066      
unitconv        = electron_charge*rydberg2ev/(1e9*bohr_radius**3)*factor

#-------------------------------------------------------------------------------

output_birch = open('birch-murnaghan',"w")

energy = [] ; volume = []

lstrain = 0.

while True:
    line = input_energy.readline().strip()
    if len(line) == 0: break
    xvolume = float(line.split()[0])
    if (lexses): 
        lstrain = xvolume
        estrain = sqrt(1.+2.*lstrain)-1
        xvolume = (1.+estrain)**dimension * vzero
    if (abs(lstrain) < abs(maxstrain)+1e-6):
        volume.append(xvolume)
        energy.append(float(line.split()[1]))

#-------------------------------------------------------------------------------

volume, energy = sortstrain(volume,energy)
nv = len(volume)
if (nv < 4): sys.exit("\n ERROR: Too few volumes ("+str(nv)+")!\n")

#-------------------------------------------------------------------------------

p = bmfit(volume,energy)

v0 = p[0][0] ; e0 = p[0][1] ; b0 = p[0][2]*unitconv ; bp =  p[0][3]

a0sc=(1*p[0][0])**(0.33333333333)
abcc=(2*p[0][0])**(0.33333333333)
afcc=(4*p[0][0])**(0.33333333333)

chi = 0 ; ebm = [] ; eee = [] ; fchi = 0
for i in range(len(volume)): 
    chi=chi+residuals(p[0],energy[i],volume[i])**2
    fchi=fchi+(energy[i]-eos(volume[i],v0,e0,b0/unitconv,bp))**2
    ebm.append(eos(volume[i],v0,e0,b0/unitconv,bp))

chi=sqrt(chi)/len(volume)

fmt='%10.4f' ; amt='%10.4f' ; emt='%9.5f' ; bmt='%8.3f' ; pmt='%16.10f' ; lmt='%10.2f'
a2t='%12.3f' ; a3t='%14.3f' ; al0='%7.4f'

string     = "     V0        B0        BP         a-sc       a-bcc      a-fcc     log(chi)"
    
if (lexses):
    a2  = b0 * bfactor
    a3  = -a2 * (dimension*(bp-2.) + 6.)  
    mls = vol2eta(v0,vzero,dimension)
    string = "        A2            A3           lagrangian strain at minimum      log(chi)"
   
print
print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
print string
if (lexsev):
    print fmt%(v0), bmt%(b0), bmt%(bp)," ",  
    print amt%(a0sc), amt%(abcc), amt%(afcc), lmt%(log10(chi))
    print >> output_birch, pmt%(v0), pmt%(b0), pmt%(bp),  
    print >> output_birch, pmt%(a0sc), pmt%(abcc), pmt%(afcc),
    print >> output_birch, pmt%(log10(chi))
else:
    if (lplanar):
       print a2t%(a2), a3t%(a3),"    ", emt%(mls), 
       print " ( alat0 =", al0%((1.+mls)*alat),") ",
    else:
       print a2t%(a2), a3t%(a3),"   ", "          ", emt%(mls), "           ",
    print  lmt%(log10(chi))
    print >> output_birch, pmt%(v0), pmt%(b0), pmt%(bp),  
    print >> output_birch, pmt%(mls), pmt%(mls), pmt%(mls),
    print >> output_birch, pmt%(log10(chi))      
print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
print

input_energy.close()
output_birch.close()

#*******************************************************************************

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])

#-------------------------------------------------------------------------------
   
xlabel = u'Volume [Bohr\u00B3]'
if (lexses): xlabel = u'Lagrangian strain'
ylabel = r'Energy [Ha]'
if (lrydberg): ylabel = r'Energy [Ry]'

#-------------------------------------------------------------------------------

xvol = numpy.linspace(volume[0],volume[-1],100)

xene = [] 
for i in range(len(xvol)):
    xene.append(eos(xvol[i],v0,e0,b0/unitconv,bp))    

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

if (lexses):
    avolume = [] ; axvol = []
    for i in range(len(volume)): avolume.append(vol2eta(volume[i],vzero,dimension))
    for i in range(len(xvol)): axvol.append(vol2eta(xvol[i],vzero,dimension))
    av0     = vol2eta(v0,vzero,dimension)
else:
    av0     = v0
    avolume = volume
    axvol   = xvol
    
plt.plot(axvol,xene,'b-',label='birch-murnaghan fit')
plt.plot(avolume,energy,'go',label='calculated')
plt.plot(av0,e0,'ro')
plt.legend(loc=9,borderaxespad=.8,numpoints=1)

ymax  = max(max(xene),max(energy))
ymin  = min(min(xene),min(energy))
dxx   = abs(max(avolume)-min(avolume))/18
dyy   = abs(ymax-ymin)/18
ax.yaxis.set_major_formatter(yfmt)
ax.set_xlim(min(avolume)-dxx,max(avolume)+dxx)
ax.set_ylim(ymin-dyy,ymax+dyy)

ax.xaxis.set_major_locator(MaxNLocator(7))

plt.savefig('PLOT.ps', orientation='portrait',format='eps')
plt.savefig('PLOT.png',orientation='portrait',format='png',dpi=dpipng)

#-------------------------------------------------------------------------------

