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

if (str(os.path.exists('input.xml'))=='False'): 
    sys.exit("ERROR: Input file input.xml not found!\n")

maximum_strain = 0.05
strain_points = input("\nEnter the number of volume values >>>> ")
strain_points = int(abs(strain_points))
if (3 > strain_points or strain_points > 99): 
    sys.exit("ERROR: Number of volume values is out of range [3-99]!\n")

#-------------------------------------------------------------------------------

input_obj = open("input.xml","r")
input_doc = etree.parse(input_obj)
input_rut = input_doc.getroot()
 
xml_scale = map(float,input_doc.xpath('/input/structure/crystal/@scale'))
if (xml_scale == []): 
    ref_scale = 1.0
else: 
    ref_scale = float(xml_scale[0])

str_stretch = input_doc.xpath('/input/structure/crystal/@stretch')
if (str_stretch ==[]): 
    xml_stretch = [1.,1.,1.]
else: xml_stretch=numpy.array(map(float,str_stretch[0].split()))

lst_basevect = input_doc.xpath('//basevect/text()')
xml_basevect = []
for ind_basevect in lst_basevect:
    xml_basevect.append(map(float,ind_basevect.split()))

axis_matrix = numpy.array(xml_basevect) 
determinant = numpy.linalg.det(axis_matrix)
volume = abs(determinant*ref_scale**3\
    *xml_stretch[0]*xml_stretch[1]*xml_stretch[2])

work_directory = 'workdir'
if (len(sys.argv) > 1): work_directory = sys.argv[1]
if (os.path.exists(work_directory)): shutil.rmtree(work_directory)
os.mkdir(work_directory)
os.chdir(work_directory)

#-------------------------------------------------------------------------------

delta=strain_points-1
convert=1
if (strain_points <= 1):
    strain_points=1
    convert=-1
    delta=1
eta_step=2*maximum_strain/delta

#-------------------------------------------------------------------------------

for i in range(0,strain_points):

#-------------------------------------------------------------------------------

    eta=i*eta_step-maximum_strain*convert

    if (i+1 < 10): strainfile = 'volume-0'+str(i+1)
    else: strainfile = 'volume-'+str(i+1)
    output_str = open(strainfile,"w")
    fmt = '%20.8f'
    vol=(1.+eta)**3*volume
    print >>output_str, (fmt%vol)
    output_str.close()

    if (abs(eta) < 0.000001): eta=0.000001

    new_scale=ref_scale*(1.+eta)

#-------------------------------------------------------------------------------

    xsc = input_doc.xpath('//crystal')
    fmt = '%16.10f'
    xsc[0].set("scale",str(fmt%new_scale))
    if (i+1 < 10): outputfile = 'input-0'+str(i+1)+'.xml'
    else: outputfile = 'input-'+str(i+1)+'.xml'
    output_obj = open(outputfile,"w")
    output_obj.write(etree.tostring(input_rut, method='xml',
                                               pretty_print=True,
                                               xml_declaration=True,
                                               encoding='UTF-8'))
    output_obj.close()

#-------------------------------------------------------------------------------

os.chdir('../')
print 

