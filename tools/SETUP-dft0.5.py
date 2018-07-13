#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   lxml  import etree
from   sys   import stdin
from   numpy import *
import subprocess
import os.path
import shutil
import numpy 
import math
import sys
import os

#-------------------------------------------------------------------------------

rcutmin   = input("\nEnter minimum value of r_cut_min [Bohr]            >>>> ")
rcutmax   = input("\nEnter maximum value of r_cut_max [Bohr]            >>>> ")
rcutsteps = input("\nEnter the number of steps in [r_cut_min,r_cut_max] >>>> ")

if (str(os.path.exists('input.xml'))=='False'): 
    sys.exit("ERROR: Input file input.xml not found!\n")

if (1 > rcutsteps or rcutsteps > 99): 
    sys.exit("ERROR: Number of r_cut values is out of range [1-99]!\n")
    
if (rcutmin >= rcutmax): rcutsteps = 1
 
#-------------------------------------------------------------------------------

input_obj = open("input.xml","r")
input_doc = etree.parse(input_obj)
input_rut = input_doc.getroot()
 
xml_cut = input_doc.xpath('/input/structure/species/dfthalfparam/@cut')
ref_cut = []
for i in range(len(xml_cut)): ref_cut.append(float(xml_cut[i])) 

#-------------------------------------------------------------------------------

work_directory = 'workdir'
if (len(sys.argv) > 1): work_directory = sys.argv[1]
if (os.path.exists(work_directory)): shutil.rmtree(work_directory)
os.mkdir(work_directory)
os.chdir(work_directory)

#-------------------------------------------------------------------------------

cstep = 1.0
if not(rcutsteps == 1): cstep = (rcutmax-rcutmin)/float(rcutsteps-1)

for i in range(0,rcutsteps):

#-------------------------------------------------------------------------------

    eta=rcutmin+i*cstep
    if (i+1 < 10): strainfile = 'strain-0'+str(i+1)
    else: strainfile = 'strain-'+str(i+1)
    output_str = open(strainfile,"w")
    fmt = '%20.8f'
    print >>output_str, (fmt%eta)
    output_str.close()
    new_cut = eta

#-------------------------------------------------------------------------------

    dft05 = input_doc.xpath('/input/structure/species/dfthalfparam')
    fmt = '%16.10f'
    dft05[-1].set("cut",str(fmt%new_cut))
    if (i+1 < 10): outputfile = 'input-0'+str(i+1)+'.xml'
    else: outputfile = 'input-'+str(i+1)+'.xml'
    output_obj = open(outputfile,"w")
    output_obj.write(etree.tostring(input_rut, method='xml',
                                               pretty_print=True,
                                               xml_declaration=False,
                                               encoding='UTF-8'))
    output_obj.close()

#-------------------------------------------------------------------------------

os.system("touch dft-0.5") 
os.chdir('../')
os.system("cp input.xml "+work_directory+"/source.xml") 
print 

