#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   lxml  import etree
from   sys   import stdin
from   numpy import *
from   pylab import *
import matplotlib        as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ptk 
import pylab             as pyl
import xml.etree.ElementTree as ET
import subprocess
import os.path
import shutil
import numpy 
import math
import sys
import xml
import os

#-------------------------------------------------------------------------------

ndir = 3

alat = [] ; fdir = [] ; qpoint = [] ; etot = []

for idir in range(ndir):
 
    ddir = raw_input('\n>>>> Directory '+str(idir+1)+' of '+str(ndir)+' : ')
    fdir.append(ddir)
    os.chdir(ddir)
    
    linput  = os.path.exists('input.xml')
    lphonon = os.path.exists('PHONON.OUT')  
    
    if (not linput): sys.exit("ERROR: "+str(ddir)+"/input.xml not found!\n")
    if (not lphonon): sys.exit("ERROR: "+str(ddir)+"/PHONON.OUT not found!\n")
  
    fileparse  = 'input.xml'    
    tree = ET.parse(fileparse)
    root = tree.getroot()

    for child in root:      
        if (child.tag == 'structure'): 
           for subchild in child:
	       if (subchild.tag == 'crystal'):
		   alat.append(subchild.get('scale')) 
      
#-------------------------------------------------------------------------------
# Read data

    os.system('CONVERT-au2invcm.py > phonon.out')
    
    f = open('phonon.out','r')
    lines = f.readlines()
    
    natom = int(lines[4].strip().split()[5])
    nqpt  = int(lines[5].strip().split()[5])
    nbnd  = natom*3
    
    ini   = 8
    eband = []
   
    for iqpt in range(nqpt):
        qpoint.append(lines[iqpt*nbnd+ini].strip())
        for ibnd in range(nbnd):
	    #print lines[iqpt*nbnd+ibnd+ini+1].strip()
            eband.append(float(lines[iqpt*nbnd+ibnd+ini+1].strip().split()[4]))
        ini=ini+2
    etot.append(eband)
    
    f.close
            
    os.chdir('../')
    
print

#-------------------------------------------------------------------------------
# Calculate Grüneisen mode parameters

gruen = []    
dalat = (float(alat[0]) - float(alat[2])) / float(alat[1]) *3.

for iqpt in range(nqpt):
    gband = []
    for ibnd in range(nbnd):
        if (abs(etot[1][iqpt*nbnd+ibnd])<0.01):
	    gband.append(0.)
        else:
            dfreq = (etot[2][iqpt*nbnd+ibnd] - etot[0][iqpt*nbnd+ibnd]) / etot[1][iqpt*nbnd+ibnd] 
            gband.append(dfreq/dalat)
    gruen.append(gband)
    
fmt = "%6.3f"
imt = "%3i"

print lines[1].strip()
print
for iqpt in range(nqpt):
    print qpoint[iqpt]
    for ibnd in  range(nbnd):
        print '  mode',imt%(ibnd),' with Grüneisen parameter ',fmt%(gruen[iqpt][ibnd])
    print
print lines[1].strip()
print

#fmt = "%18.14f"

#outfile = open('GRDISP.OUT',"w")

#for ibnd in range(nbnd):
#    for i in range(nkpt):
#        print >>outfile, fmt%(kline[i]), fmt%(gruen[ibnd][i])
#    print >>outfile

#outfile.close()

sys.exit()    

#-------------------------------------------------------------------------------

