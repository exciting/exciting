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
              input("\nEnter minimum interlayer distance dmin [Bohr] >>>> ")
maximum_displ = \
              input("\nEnter maximum interlayer distance dmax [Bohr] >>>> ")
displ_points  = \
              input("\nEnter the number of distances in [dmin,dmax]  >>>> ")
infty_displ_s = \
              raw_input("\nEnter interlayer distance at infinity  [Bohr] >>>> ")

linfty = True                         
if ((str(infty_displ_s) == "")): 
    linfty = False
else:
    if ((float(infty_displ_s) <= maximum_displ )): linfty = False
if (linfty): infty_displ = float(infty_displ_s)
             
displ_points = int(abs(displ_points))
if (0 > displ_points or displ_points > 99): 
    sys.exit("ERROR: Number of displacements is out of range [1-99]!\n")

#-------------------------------------------------------------------------------

input_obj = open("input.xml","r")
input_doc = etree.parse(input_obj)
input_rut = input_doc.getroot()
  
xml_scale = map(float,input_doc.xpath('/input/structure/crystal/@scale'))
if (xml_scale == []): ref_scale = 1.0
else: ref_scale = float(xml_scale[0])

xml_stretch = input_doc.xpath('/input/structure/crystal/@stretch')
if (xml_stretch ==[]): ref_stretch = [1.,1.,1.]
else: ref_stretch = numpy.array(map(float,xml_stretch[0].split()))

xml_basevect = input_doc.xpath('//basevect/text()')
ref_basevect = []
for ind_basevect in xml_basevect: ref_basevect.append(map(float,ind_basevect.split()))
mat_basevect = numpy.array(ref_basevect) 
new_basevect = mat_basevect
 
xml_atom = input_doc.xpath('/input/structure/species/atom/@coord')
ref_atom = []
for ind_atom in xml_atom: ref_atom.append(map(float,ind_atom.split()))
mat_atom = numpy.array(ref_atom) 

#print ref_scale
#print ref_stretch
#print mat_basevect
#print new_basevect
#print mat_atom

work_directory = 'workdir'
if (len(sys.argv) > 1): work_directory = sys.argv[1]
if (os.path.exists(work_directory)): shutil.rmtree(work_directory)
os.mkdir(work_directory)
os.chdir(work_directory)

#-------------------------------------------------------------------------------

delta      = displ_points-1
displ_step = float(maximum_displ-minimum_displ)/delta

#-------------------------------------------------------------------------------

for i in range(0,displ_points):

#-------------------------------------------------------------------------------

    displ = minimum_displ + i*displ_step
    
    if (i+1 < 10): displfile = 'strain-0'+str(i+1)
    else: displfile = 'strain-'+str(i+1)
    output_dsp = open(displfile,"w")
    fmt = '%11.8f'
    print >>output_dsp, (fmt%displ)
    output_dsp.close()
    
    new_basevect[2,2] = 2.*displ/ref_scale/ref_stretch[2]
    
#-------------------------------------------------------------------------------

    xbv = input_doc.xpath('//crystal/basevect')
    fmt = '%22.16f'
    for j in range(3):
        xbv[j].text = str(fmt%new_basevect[j,0])\
                    + str(fmt%new_basevect[j,1])\
                    + str(fmt%new_basevect[j,2]) +" "
    if (i+1 < 10): outputfile = 'input-0'+str(i+1)+'.xml'
    else: outputfile = 'input-'+str(i+1)+'.xml'
    output_obj = open(outputfile,"w")
    output_obj.write(etree.tostring(input_rut, method='xml',
                                               pretty_print=True,
                                               xml_declaration=True,
                                               encoding='UTF-8'))
    output_obj.close()

#-------------------------------------------------------------------------------

if (linfty):
    displ = infty_displ
    displfile = 'strain-oo'
    output_dsp = open(displfile,"w")
    fmt = '%11.8f'
    print >>output_dsp, (fmt%displ)
    output_dsp.close()
    
    new_basevect[2,2] = 2.*displ/ref_scale/ref_stretch[2]
    
    xbv = input_doc.xpath('//crystal/basevect')
    fmt = '%22.16f'
    for j in range(3):
        xbv[j].text = str(fmt%new_basevect[j,0])\
                    + str(fmt%new_basevect[j,1])\
                    + str(fmt%new_basevect[j,2]) +" "
    outputfile = 'input-oo.xml'
    output_obj = open(outputfile,"w")
    output_obj.write(etree.tostring(input_rut, method='xml',
                                               pretty_print=True,
                                               xml_declaration=True,
                                               encoding='UTF-8'))
    output_obj.close()

#-------------------------------------------------------------------------------

output_label = open("xlabel","w")
print >>output_label, "Interlayer distance [Bohr]"
output_label.close()

#-------------------------------------------------------------------------------

os.chdir('../')
print 

