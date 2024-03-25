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

def mantissa(x):
    logx = log10(abs(x)) 
    espo = int(logx)
    newx = 10**(logx - espo)
    while newx < 1:
        newx = newx*10 ; espo = espo-1
    return newx, espo 

#-------------------------------------------------------------------------------

def purifica(x,y,ymin,ymax):
    checkp = 0
    j = -1
    for i in range(len(y)):
        if (ymax < y[j] or ymin > y[j]): 
            del y[j]
            del x[j]     
        else:
            checkp = 1
        if (checkp != 0): break
    return

#-------------------------------------------------------------------------------

def mymantissa(x):
    a,b = mantissa(x) ; mant = 10.**b
    fmt4  = '%6.4f' ; f4 = fmt4%(abs(x)/mant)
    ss = '10^'+str(b)
    if (b == 0): ss = ''
    if (b == 1): ss = '$10\,$'
    sx = u'\u2013'+f4+'*'+ss
    if (x > 0): sx = u'+'+f4+'*'+ss
    return sx   
    
#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

#if (narg<2): 
#    print "\nIncorrect number of arguments. **Usage**:\n\n",
#    print "PLOT-checkderiv.py YMIN YMAX\n"
#    sys.exit()

ylimits = []
if (len(sys.argv) > 1): ylimits.append(float(sys.argv[1]))
if (len(sys.argv) > 2): ylimits.append(float(sys.argv[2]))

#-------------------------------------------------------------------------------

unit   = r'GPa]'
xlabel = r'Maximum lagrangian strain'

if (os.path.exists('energy-vs-strain')): 
    unit   = r'GPa]'
    xlabel = r'Maximum lagrangian strain'
    if (os.path.exists('planar')): unit = r'N/m]'

if (os.path.exists('energy-vs-displacement')): 
    unit   = u'cm$^-\!$\u00B9]'
    xlabel = r'Maximum displacement $u$ [alat]' 

inpf='check-energy-derivatives'

if (str(os.path.exists('order-of-derivative'))=='False'): 
    sys.exit("\nERROR: file \"order-of-derivative\" not found!\n")

startorder = 0 
if (os.path.exists('startorder')): 
    startorder = int(open('startorder',"r").readline().strip().split()[0])

order  = int(open('order-of-derivative',"r").readline().strip().split()[0])
ylabel = r'Derivative of order $'+str(order)+'$'
lg     = order+startorder 

if (str(os.path.exists(inpf))=='False'): 
    sys.exit("\nERROR: file "+inpf+" not found!\n")
    
#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])
pdefault = shell_value('PDEFAULT',ev_list,"")[1] 
writeminmax = shell_value('WRITEMINMAX',ev_list,"")[0]

#-------------------------------------------------------------------------------
   
input_file = open(inpf,"r")

x1 = [] ; x2 = [] ; x3 = [] 
y1 = [] ; y2 = [] ; y3 = []

while True:
    line = input_file.readline()
    line = line.strip()
    if len(line) == 0: break
    if (len(line.split()) > 1): 
        y1.append(float(line.split()[1])) 
        x1.append(float(line.split()[0]))
    if (len(line.split()) > 3): 
        y2.append(float(line.split()[3]))
        x2.append(float(line.split()[0]))
    if (len(line.split()) > 5): 
        y3.append(float(line.split()[5]))
        x3.append(float(line.split()[0]))

#-------------------------------------------------------------------------------
# manipulate data for a better plot

if (len(y1) > 0): rmin  = y1[0] 
if (len(y2) > 0): rmin  = y2[0] 
if (len(y3) > 0): rmin  = y3[0] 

newrmin,espo = mantissa(rmin) ; mant = 10.**espo

fmt4  = '%6.4f' ; ffour = fmt4%(abs(rmin)/mant)

ssexp = '$10^{'+str(espo)+'}$'
if (espo == 0): ssexp = ''
if (espo == 1): ssexp = '$10\,$'
luni  = ffour+' ['+ssexp+unit
srmin = u'\u2013 '+luni
if (rmin > 0): srmin = u'+ '+luni

for i in range(len(y1)): y1[i]=(y1[i]-rmin)/mant
for i in range(len(y2)): y2[i]=(y2[i]-rmin)/mant
for i in range(len(y3)): y3[i]=(y3[i]-rmin)/mant

ymin  = min(y1+y2+y3) ; ymax  = max(y1+y2+y3)

dyy   = abs(ymax-ymin)/18 ; ymin  = ymin-dyy ; ymax = ymax+dyy

if (len(ylimits) == 1): 
    ymin = (float(ylimits[0])-rmin)/mant
if (len(ylimits) > 1): 
    ymin = (float(ylimits[0])-rmin)/mant ; ymax = (float(ylimits[1])-rmin)/mant 

purifica(x1,y1,ymin,ymax)
purifica(x2,y2,ymin,ymax)
purifica(x3,y3,ymin,ymax)

xmin  = min(x1+x2+x3) ; xmax  = max(x1+x2+x3)
dxx   = abs(xmax-xmin)/18 ; xmin  = xmin-dxx ; xmax = xmax+dxx

#-------------------------------------------------------------------------------
# set defauls parameters for the plot

fontlabel=20
fonttick=16
fonttext=18
fontlimits=12

params = {'ytick.minor.size': 6,
          'xtick.major.pad': 8,
          'ytick.major.pad': 4,
          'patch.linewidth': 2.,
          'axes.linewidth': 2.,
          'lines.linewidth': 1.8,
          'lines.markersize': 8.0,
          'mathtext.fontset' : 'stixsans',
          'axes.formatter.limits': (-8, 8)}

plt.rcParams.update(params)
plt.subplots_adjust(left=0.20, right=0.78,
                    bottom=0.18, top=0.88,
                    wspace=None, hspace=None)
                    
yfmt = ptk.ScalarFormatter(useOffset=True,useMathText=True)

fig  = matplotlib.pyplot.figure(1, figsize=(8,5.5)) 

ax   = fig.add_subplot(111)

ax.text(0.5,-0.13,xlabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=0)
#ax.text(-0.19,0.5,ylabel,size=fontlabel,
#        transform=ax.transAxes,ha='center',va='center',rotation=90)

ax.text(0.0,1.05,srmin,size=fonttext,color='#00008B',
        transform=ax.transAxes,ha='left',va='center',rotation=0)
        
if (writeminmax=="1"): 
    ax.text(0.62,1.033,"m = "+mymantissa(ymin*mant+rmin),size=fontlimits,color='#00008B',
        transform=ax.transAxes,ha='left',va='center',rotation=0)
 
    ax.text(0.62,1.080,"M = "+mymantissa(ymax*mant+rmin),size=fontlimits,color='#00008B',
        transform=ax.transAxes,ha='left',va='center',rotation=0)
 
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(6)
    line.set_markeredgewidth(2)      
       
plt.xticks(size=fonttick)
plt.yticks(size=fonttick)
pyl.grid(True)
plt.ylabel(ylabel,size=fontlabel)

if (len(y1) > 0): plt.plot(x1,y1,'ro-',label='n='+str(lg))
if (len(y2) > 0): plt.plot(x2,y2,'bo-',label='n='+str(lg+2))
if (len(y3) > 0): plt.plot(x3,y3,'go-',label='n='+str(lg+4))
plt.legend(bbox_to_anchor=(1.03, 1), loc=2, borderaxespad=0.)

ax.yaxis.set_major_formatter(yfmt)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

ax.xaxis.set_major_locator(MaxNLocator(7))

ax.set_axisbelow(True) 

plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------





