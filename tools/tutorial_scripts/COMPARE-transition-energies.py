#!/usr/bin/python2
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

#-------------------------------------------------------------------------------
# Create the list of input directories 

test=True
dirs=[]

print("\n################################################\n")
print(" Enter the names of the working directories and ")
print(' enter "(Q)uit" or "(q)uit" to terminate input.\n')

while test:
	dir=raw_input(" Directory name ==> ")
	if dir[0:1]=="q" or dir[0:1]=="Q" or dir=="Quit" or dir=="quit": test=False
	dirs.append(dir)
	

GGgap = []
GXgap = []
for j in range(len(dirs)-1):
	dir=dirs[j]
	infile=root+"/"+dir+"/EIGVAL.OUT"
	bashCommand = "grep \"nstsv\" "+infile+" |awk '{print $1}'"
	#bashCommand = "grep -A56  \"0.000000000       0.000000000       0.000000000\" "+infile
	stream=os.popen(bashCommand)
	nst=stream.read()
	bashCommand = "grep -A"+ str(int(nst)+1) + "  \"0.000000000       0.000000000       0.000000000\" " \
                      +infile + "| grep \"2.000000000\"| awk 'BEGIN{a=0}{if ($2>0+a) a=$2} END{print a*27.211}' " 
	GammaVBM=os.popen(bashCommand).read()
	bashCommand = "grep -A"+ str(int(nst)+1) + "  \"0.000000000       0.000000000       0.000000000\" " \
                      +infile +"| tail -n"+ str(int(nst)-1)+ "| grep \"0.000000000\"| awk 'BEGIN{a=1000}\
                      {if ($2<0+a) a=$2} END{print a*27.211}' "
	GammaCBm=os.popen(bashCommand).read()
	bashCommand = "grep \"0.5000000000      0.5000000000       0.000000000\" " +infile 
	stream=os.popen(bashCommand).read()  
	if stream != '':
		bashCommand = "grep -A"+ str(int(nst)+1) + "  \"0.5000000000      0.5000000000       0.000000000\" "\
                              +infile +"| tail -n"+ str(int(nst)-1)+ "| grep \"0.000000000\"| awk 'BEGIN{a=1000}\
                              {if ($2<0+a) a=$2} END{print a*27.211}' "
		XCBm=os.popen(bashCommand).read()
	else:
		bashCommand = "grep -A"+ str(int(nst)+1) + "  \"0.000000000      0.5000000000      0.5000000000\" " \
                              +infile +"| tail -n"+ str(int(nst)-1)+ "| grep \"0.000000000\"| awk 'BEGIN{a=1000}\
                              {if ($2<0+a) a=$2} END{print a*27.211}' "
		XCBm=os.popen(bashCommand).read()
	GGgap.append(float(GammaCBm)-float(GammaVBM))
	GXgap.append(float(XCBm)-float(GammaVBM))
print("\n------------------------------------------------\n")
print " Transition energies in eV:\n"
fmtL = "       ".join(["{:.6s}"]*len(GGgap))
print "                   ",   fmtL.format(*dirs[:len(dirs)-1] )
fmtL = "     ".join(["{:.3f}"]*len(GGgap))
print " Gamma -> Gamma:   ", fmtL.format(*GGgap[:])   
print " Gamma -> X    :   ", fmtL.format(*GXgap[:])   
print("\n################################################")
	
