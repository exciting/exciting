#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   subprocess import *
from   lxml       import etree as ET
from   scipy.optimize.optimize import fmin
from   scipy.optimize.optimize import fmin_powell
import numpy as np
import time
import sys
import os
  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def include_relax(ifile,gfile):
  
    lfile = os.path.exists(gfile)
    if (not(lfile)):  return 
  
    itree = ET.parse(ifile)
    gtree = ET.parse(gfile)
    root  = itree.getroot()
    
    ixml_atom    = itree.xpath('//species/atom')
    gxml_species = gtree.xpath('//species')
    
    amt   = '%15.9f'
    katom = 0
    for i in range(len(gxml_species)):
        gxml_atom = gxml_species[i].findall('atom')
        for j in range(len(gxml_atom)):  
            atom = np.array(map(float,gxml_atom[j].get("coord").split()))
            line = (amt%atom[0]+amt%atom[1]+amt%atom[2])
            ixml_atom[katom].set("coord",line)
            katom = katom+1
            
    OUTOBJ = open('opt-rel.xml', 'w')
    OUTOBJ.write(ET.tostring(root, method         ='xml',
                                   pretty_print   =True ,
                                   xml_declaration=False ,
                                   encoding       ='UTF-8'))
    OUTOBJ.close()
    
    return 

#START==========================================================================          

narg  = len(sys.argv)-1

ifile = "input.xml"
gfile = "geometry_opt.xml"

if (narg>0): ifile = str(sys.argv[1])
if (narg>1): gfile = str(sys.argv[2])

lifile = os.path.exists(ifile)
lgfile = os.path.exists(gfile)

if (not(lifile)): sys.exit('\n ERROR: There is NO '+ifile+' file !?!?!? \n')

os.system("cp "+ifile+" "+ifile[0:len(ifile)-4]+"_rel.xml")

if (lgfile): 
    include_relax(ifile, gfile)
    os.system("mv opt-rel.xml "+ifile[0:len(ifile)-4]+"_rel.xml")
   
#END============================================================================          
         

