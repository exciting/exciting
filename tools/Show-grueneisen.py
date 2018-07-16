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

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

#-------------------------------------------------------------------------------

ndir = 3

alat = [] ; fdir = [] ; eband = []

for idir in range(ndir):
 
    ddir = raw_input('\n>>>> Directory '+str(idir+1)+' of '+str(ndir)+' : ')
    fdir.append(ddir)
    os.chdir(ddir)
    
    linput  = os.path.exists('source.xml')
    
    if (not linput): sys.exit("ERROR: "+str(ddir)+"/input.xml not found!\n")

    fileparse  = 'source.xml'    
    tree = ET.parse(fileparse)
    root = tree.getroot()

    for child in root:      
        if (child.tag == 'structure'): 
           for subchild in child:
	       if (subchild.tag == 'crystal'):
		   alat.append(subchild.get('scale')) 
      
#-------------------------------------------------------------------------------
# Read data
    
    f = open('phonon.out','r')
    lines = f.readlines()
    eband.append(float(lines[0].strip().split()[0]))
    f.close
    os.chdir('../')

#-------------------------------------------------------------------------------
# Calculate Grüneisen mode parameters

gruen = []    
dalat = (float(alat[0]) - float(alat[2])) / float(alat[1]) *3.

if (abs(eband[1])<0.01):
    gruen = 0.
else:
    dfreq = (eband[2] - eband[0]) / eband[1] 
    gruen = dfreq / dalat
    
fmt = "%6.3f"
imt = "%3i"

print
print "#############################################"
print
print '  mode Grüneisen parameter :          ',fmt%(gruen)
print
print "#############################################"
print

sys.exit()    

#-------------------------------------------------------------------------------

