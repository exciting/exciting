#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   lxml  import etree
from   sys   import stdin
from   numpy import *
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

def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1
ndir = 3
if (narg>0): ndir = int(sys.argv[1])
 
#-------------------------------------------------------------------------------

alat = [] ; fdir = []

for idir in range(ndir):
 
    ddir = raw_input('\n>>>> Directory '+str(idir+1)+' of '+str(ndir)+' : ')
    fdir.append(ddir)
    os.chdir(ddir)
    
    linput     = os.path.exists('input.xml')
    linput_ini = os.path.exists('input_ini.xml')
    
    if (not linput): sys.exit("ERROR: "+str(ddir)+"/input.xml not found!\n")
   
    fileparse  = 'input_ini.xml'    
    if (not linput_ini): fileparse  = 'input.xml' 
          
    tree = ET.parse(fileparse)
    root = tree.getroot()
    os.rename(fileparse,'input_ini.xml')

    for child in root:
      
        if (child.tag == 'structure'): 
           for subchild in child:
	       if (subchild.tag == 'crystal'):
		   alat.append(subchild.get('scale')) 
      
        if (child.tag == 'groundstate'): 
           child.set('do', 'skip')
        
        if (child.tag == 'phonons'): 
           child.set('do', 'skip')
           
           for subchild in child: child.remove(subchild)
        
           phonondos = ET.SubElement(child, 'phonondos')
           phonondos.set('ngrdos', '100')
           phonondos.set('nwdos', '500')
           phonondos.set('ntemp', '200')
           phonondos.set('nsmdos', '2') 
            
           phonondispplot = ET.SubElement(child, 'phonondispplot')
           plot1d = ET.SubElement(phonondispplot, 'plot1d')
           path = ET.SubElement(plot1d, 'path')
           path.set('steps', '200')
           point = ET.SubElement(path, 'point')
           point.set('coord', '1.0   0.0   0.0')
           point = ET.SubElement(path, 'point')
           point.set('coord', '0.5   0.5   0.0')
           point = ET.SubElement(path, 'point')
           point.set('coord', '0.0   0.0   0.0')
           point = ET.SubElement(path, 'point')
           point.set('coord', '0.5   0.0   0.0')
 
    indent(root)    
    tree.write('input.xml')  
    
    os.system('touch phrun')
    
    os.chdir('../')

#-------------------------------------------------------------------------------
 
info = open('INFO-thermox',"w")

print >>info, "\nTotal number of directories = ",ndir
print >>info

for idir in range(ndir):
    print >>info, "Cubic lattice parameter = ",
    print >>info, str(alat[idir])+" [a.u.]   # Directory = ",
    print >>info, str('"')+str(fdir[idir])+str('"')
print >>info

info.close()

print

#-------------------------------------------------------------------------------

