#!/usr/bin/python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# -3c               3d surface plot
# -3w               3d wire plot
# -3c               3d contour plot
# -c                2d contour plot
# -file FILENAME    Input file is FILENAME in local directory
#_______________________________________________________________________________

import matplotlib
import matplotlib.pyplot as plt
from   matplotlib.colors import LinearSegmentedColormap
from   mpl_toolkits.mplot3d.axes3d import Axes3D
from   matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.cm as cm 
import parse_plot2d
import numpy as np
import os
import sys
from   pylab import *

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

showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])

#-------------------------------------------------------------------------------

args  = sys.argv[1:]

fname = "VCL2d2D.XML"

if "-file" in args:
    ind = args.index("-file")
    fname = str(args[ind+1]) 

if (str(os.path.exists(current+'/'+fname))=='False'): 
    sys.exit("\n ERROR: Input file "+current+"/"+fname+" not found!\n")
    
rx, ry, basevect, rz = parse_plot2d.parse_plot2d("input.xml", fname)

x = [] ; y = [] ; z = []

xgrid = len(rx)
ygrid = len(rx[1])

#print xgrid, ygrid

for i in range(xgrid):
    x.append(i*1.0/(xgrid-1))
    y.append(i*1.0/(ygrid-1))
    
Y, X = np.meshgrid(x, y)

Z=[]

zmin=1.e30
zmax=0.

for i in range(xgrid):
    d = []
    for j in range(ygrid):
        d.append(rz[i][j])
        if (abs(rz[i][j]) > zmax): zmax = abs(rz[i][j]) 
        if (abs(rz[i][j]) < zmin): zmin = abs(rz[i][j])         
    Z.append(d)  

fig = plt.figure()
ax = fig.gca(projection='3d')

if   "-c" in args:  
    fig, ax = plt.subplots()
    p = ax.pcolor(X, Y, Z, cmap=cm.coolwarm, vmin=zmin, vmax=zmax)
    cb = fig.colorbar(p, ax=ax)
elif "-3c" in args:
    ppp = ax.contour(X, Y, Z, zdir='z', offset=0.00, cmap=cm.jet)
    ax.contour(X, Y, Z, zdir='x', offset=0.00)#, cmap=cm.jet)
    ax.contour(X, Y, Z, zdir='y', offset=1.00)#, cmap=cm.jet)
    fig.colorbar(ppp, shrink=0.5, aspect=8)
    ax.zaxis.set_major_locator(LinearLocator(8))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))   
elif "-3w" in args:
    ppp = ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)
    ax.zaxis.set_major_locator(LinearLocator(8))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
elif "-3s" in args:
    ppp = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
          linewidth=0, antialiased=False)
    fig.colorbar(ppp, shrink=0.5, aspect=8)
else:
    ppp = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
               linewidth=0, antialiased=False)
    fig.colorbar(ppp, shrink=0.5, aspect=8)
    ax.zaxis.set_major_locator(LinearLocator(8))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------











