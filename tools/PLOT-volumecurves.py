#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   lxml import etree
from   numpy import *
import os
import copy
import matplotlib as mpl
import numpy
import sys

import matplotlib.transforms as ptf
import matplotlib.ticker     as ptk
import matplotlib.pyplot     as plt
import pylab                 as pyl

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

#-------------------------------------------------------------------------------

root=os.getcwd()
order_of_fit=2

#-------------------------------------------------------------------------------
# Create the list of input directories 

test=True
dirs=[]

print "\n################################################\n"
print " Enter the names of the working directories. \n"
print ' Enter "(Q)uit" or "(q)uit" to terminate input\n'
print "------------------------------------------------\n"

while test:
	dir=raw_input(" Directory name ==> ")
	if dir[0:1]=="q" or dir[0:1]=="Q" or dir=="Quit" or dir=="quit": test=False
	dirs.append(dir)
	
print "\n################################################\n"

#-------------------------------------------------------------------------------
# settings for the plot 

figcolor = 'white'
dpi = 300
fig = plt.figure(figsize=(14.5,10),dpi=dpi)
fig.figurePatch.set_edgecolor(figcolor)
fig.figurePatch.set_facecolor(figcolor)

mpl.rcParams['axes.linewidth' ] = 4.0 
mpl.rcParams['grid.linewidth' ] = 1.5
mpl.rcParams['xtick.labelsize'] = 30
mpl.rcParams['ytick.labelsize'] = 30
mpl.rcParams['axes.edgecolor' ] = 'black'
mpl.rcParams['axes.labelsize' ] = 45     
mpl.rcParams['axes.labelcolor'] = 'black'
mpl.rcParams['axes.axisbelow' ] = 'True'   
mpl.rcParams['legend.fontsize'] = 30
plt.rcParams['xtick.major.pad'] = 20
plt.rcParams['ytick.major.pad'] = 10

colors=['b','r','g','y','k']

ax1 = fig.add_axes([0.2,0.18,0.75,0.76])

ax1.xaxis.set_label_position('bottom')
ax1.set_ylabel('Energy - E$\mathregular{_{min}}$ [Ha]', labelpad=19)
ax1.set_xlabel('Volume [Bohr$\mathregular{^{3}}$]', labelpad=13)
pyl.grid(True)

#-------------------------------------------------------------------------------
# Read inputs 

for j in range(len(dirs)-1):
	dir=dirs[j]
	infile=root+"/"+dir+"/input-01.xml"
	liste=[]
	try:
		tree = etree.parse(infile)
		groundstate=tree.xpath('/input/groundstate')
		#
		# get the name of the functional ---------------------------------------
		#
		try:
			libxc=tree.xpath('/input/groundstate/libxc')
			try:
				functional=libxc[0].attrib['correlation']+"+"\
                          +libxc[0].attrib['exchange']
			except:
				functional=libxc[0].attrib['xc']
		except:
			try:
				functional=groundstate[0].attrib['xctype']
			except:
				functional="GGA_PBE_SOL"
		#		
		# read values for the plot ---------------------------------------------
		#
		input_energy = open(root+"/"+dir+"/energy-vs-volume","r")
		lines=input_energy.readlines()
		input_energy.close()
		for i in range(len(lines)):
			energy=(float(lines[i].split()[1]))
			strain=(float(lines[i].split()[0]))
			liste.append([strain,energy])
		liste=sorted(liste)
		#
		# prepare values for the plot ------------------------------------------
		#
		energy=[liste[i][1] for i in range(len(liste))]
		emin=min(energy)
		for i in range(len(energy)): energy[i]=energy[i]-emin 
		strain=[liste[i][0] for i in range(len(liste))]
		ax1.plot(strain,energy,color=colors[j],label=functional,\
                         marker='o',markersize=12,linewidth=3.0,)
	except:
		print "no valid input file in: ",dir 

#-------------------------------------------------------------------------------
# Plot

for line in ax1.get_xticklines() + ax1.get_yticklines():
    line.set_markersize(8)
    line.set_markeredgewidth(3)

leg=ax1.legend(loc="upper center",borderaxespad=0.5,numpoints=1)
leg.get_frame().set_linewidth(4.0)
leg.get_frame().set_edgecolor("grey")
leg.draw_frame(True)

plt.ylim(ymin=-0.001)
xmin, xmax = plt.xlim()
plt.hlines(0.0,xmin,xmax,linewidth=3.0,linestyles="dashed")

plt.savefig('PLOT.eps', orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=300)

sys.exit()    

#-------------------------------------------------------------------------------

