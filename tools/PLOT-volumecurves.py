#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from lxml import etree
from numpy import *
import os
import copy
import matplotlib as mpl
import numpy
import sys
#import os

import matplotlib.transforms as ptf
import matplotlib.ticker     as ptk
import matplotlib.pyplot     as plt
import pylab                 as pyl

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

root=os.getcwd()
order_of_fit=2

## Create the list of input directories ##
test=True
dirs=[]
print "\n################################################\n"
print "    Enter list of working directories "
print '    Entering "Quit" will terminate the input\n'
print "################################################\n"

while test:
	dir=raw_input("directory name: ")
	if len(dir)==0 or dir=="Quit" or dir=="quit":
		test=False
	dirs.append(dir)

## settings for the plot ##
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
mpl.rcParams['axes.labelsize'] = '25'     # fontsize of the x any y labels
mpl.rcParams['axes.labelcolor'] = 'black'
mpl.rcParams['axes.axisbelow'] = 'True'   # whether axis gridlines and ticks are below
                                          # the axes elements (lines, text, etc)
mpl.rcParams['legend.fontsize'] = '25'
plt.rcParams['xtick.major.pad'] = '15'
plt.rcParams['ytick.major.pad'] = '10'
colors=['b','r','g','y','k']
fontlabel=30
plt.subplots_adjust(left=0.16, right=0.97,
                    bottom=0.15, top=0.92,
                    wspace=None, hspace=None)
ax   = fig.add_subplot(111)
xlabel=u'Volume [Bohr\u00B3]'
ylabel=r'Energy [Ha]'
ax.text(0.5,-0.12,xlabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=0)
ax.text(-0.145,0.5,ylabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=90)
pyl.grid(True)

## Read inputs ##
for j in range(len(dirs)-1):
	dir=dirs[j]
	infile=root+"/"+dir+"/input-01.xml"
	liste=[]
	try:
		tree = etree.parse(infile)
		groundstate=tree.xpath('/input/groundstate')
		## get the name of the functional ##
		try:
			libxc=tree.xpath('/input/groundstate/libxc')
			try:
				functional=libxc[0].attrib['correlation']+"+"+libxc[0].attrib['exchange']
			except:
				functional=libxc[0].attrib['xc']
		except:
			try:
				functional=groundstate[0].attrib['xctype']
			except:
				functional="GGA_PBE_SOL"
		## read values for the plot ##
		input_energy = open(root+"/"+dir+"/energy-vs-volume","r")
		lines=input_energy.readlines()
		input_energy.close()
		for i in range(len(lines)):
			energy=(float(lines[i].split()[1]))
			strain=(float(lines[i].split()[0]))
			liste.append([strain,energy])
		liste=sorted(liste)
		## prepare values for the plot ##
		energy=[liste[i][1] for i in range(len(liste))]
                strain=[liste[i][0] for i in range(len(liste))]
                fitr = numpy.polyfit(strain,energy,order_of_fit)
                curv = numpy.poly1d(fitr)
                vmin = numpy.roots(numpy.polyder(fitr))
                xvol = numpy.linspace(strain[0],strain[-1],100)
                plt.plot(xvol,curv(xvol)-curv(vmin[0]),label=functional,color=colors[j],lw=4.0)
	except:
		print "no valid input file in: ",dir 

## plot ##
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(8)
    line.set_markeredgewidth(3)
plt.legend(loc='upper center')
#plt.legend(bbox_to_anchor=(0.85,0.95))
plt.savefig('XC.eps',  orientation='portrait',format='eps')
plt.savefig('XC.png', orientation='portrait',format='png',dpi=300)
plt.show()



