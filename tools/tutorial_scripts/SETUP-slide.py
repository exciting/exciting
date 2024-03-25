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

minimum_displ = \
              input("\nEnter minimum coordinate umin >>>> ")
maximum_displ = \
              input("\nEnter maximum coordinate umax >>>> ")
displ_points  = \
              input("\nEnter the number of displacements in [umin,umax] >>>> ")
zeta          = \
              input("\nEnter internal displacement zeta >>>> ")

displ_points = int(abs(displ_points))
if (0 > displ_points or displ_points > 99): 
    sys.exit("ERROR: Number of displacements is out of range [1-99]!\n")

#-------------------------------------------------------------------------------

input_obj = open("input.xml","r")
input_doc = etree.parse(input_obj)
input_rut = input_doc.getroot()
 
lst_atom = input_doc.xpath('/input/structure/species/atom/@coord')
xml_atom = []
for ind_atom in lst_atom:
    xml_atom.append(map(float,ind_atom.split()))
otm = numpy.array(xml_atom) 

work_directory = 'workdir'
if (len(sys.argv) > 1): work_directory = sys.argv[1]
if (os.path.exists(work_directory)): shutil.rmtree(work_directory)
os.mkdir(work_directory)
os.chdir(work_directory)

#-------------------------------------------------------------------------------

one3 = 0.333333333333333333333333
two3 = 0.666666666666666666666667

delta=displ_points-1

displ_step=(maximum_displ-minimum_displ)/delta

#-------------------------------------------------------------------------------

for i in range(0,displ_points):

#-------------------------------------------------------------------------------

    displ=minimum_displ+i*displ_step

    if (i+1 < 10): displfile = 'strain-0'+str(i+1)
    else: displfile = 'strain-'+str(i+1)
    output_dsp = open(displfile,"w")
    fmt = '%11.8f'
    print >>output_dsp, (fmt%displ)
    output_dsp.close()

    ntm = mat([[ 0.0,             0.0,             0.0], 
               [ one3+zeta,       one3+zeta,       0.0],  
               [ displ,           displ,           0.5],
               [ one3+displ+zeta, one3+displ+zeta, 0.5]])
    
#-------------------------------------------------------------------------------

    xtm=input_doc.xpath('//species/atom')
    fmt='%16.10f'
    for j in range(4): 
        xtm[j].\
              set("coord",str(fmt%ntm[j,0])+str(fmt%ntm[j,1])+str(fmt%ntm[j,2]))
    if (i+1 < 10): outputfile = 'input-0'+str(i+1)+'.xml'
    else: outputfile = 'input-'+str(i+1)+'.xml'
    output_obj = open(outputfile,"w")
    output_obj.write(etree.tostring(input_rut, method='xml',
                                               pretty_print=True,
                                               xml_declaration=False,
                                               encoding='UTF-8'))
    output_obj.close()

#-------------------------------------------------------------------------------

os.chdir('../')
print 

