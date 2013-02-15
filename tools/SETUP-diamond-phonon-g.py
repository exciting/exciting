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

maximum_displ = \
           input("\nEnter maximum displacement umax [u/(lat. parameter)] >>>> ")
if (1 < maximum_displ or maximum_displ < 0): 
    sys.exit("ERROR: Maximum displacement is out of range [0-1]!\n")
displ_points = \
              input("\nEnter the number of displacements in [-umax,umax] >>>> ")
displ_points = int(abs(displ_points))
if (3 > displ_points or displ_points > 99): 
    sys.exit("ERROR: Number of displacements is out of range [3-99]!\n")

#-------------------------------------------------------------------------------

input_obj = open("input.xml","r")
input_doc = etree.parse(input_obj)
input_rut = input_doc.getroot()
 
xml_scale = map(float,input_doc.xpath('/input/structure/crystal/@scale'))
if (xml_scale == []): 
    sys.exit("ERROR: There is NO scale attribute in input.xml!\n")

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
volume = abs(determinant*xml_scale[0]**3\
    *xml_stretch[0]*xml_stretch[1]*xml_stretch[2])

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

output_info = open('INFO-diamond-phonon',"w")
print >>output_info, "\nEquilibrium lattice parameter (alat) = ",\
                                             ('%11.8f'%xml_scale[0]), "[a.u.]",\
                     "\nMaximum displacement (u,u,u);      u = ",\
                                                       maximum_displ, "[alat]",\
                     "\nNumber of displacements              = ",\
                                                                  displ_points,\
                     "\nVolume of equilibrium unit cell      = ",\
                                                         volume, "[a.u]^3", "\n"
output_info.close()

#-------------------------------------------------------------------------------

delta=displ_points-1
convert=1
if (displ_points <= 1):
    displ_points=1
    convert=-1
    delta=1

displ_step=2*maximum_displ/delta

#-------------------------------------------------------------------------------

for i in range(0,displ_points):

#-------------------------------------------------------------------------------

    displ=i*displ_step-maximum_displ*convert

    if (i+1 < 10): displfile = 'displ-0'+str(i+1)
    else: displfile = 'displ-'+str(i+1)
    output_dsp = open(displfile,"w")
    fmt = '%11.8f'
    print >>output_dsp, (fmt%displ)
    output_dsp.close()

    if (abs(displ) < 0.000001): displ=0.000001
    uuu = mat([[ 0.0, 0.0, 0.0], [ displ, displ, displ]])
    ntm = otm + uuu

#-------------------------------------------------------------------------------

    xtm=input_doc.xpath('//species/atom')
    fmt='%16.10f'
    for j in range(2): 
        xtm[j].\
              set("coord",str(fmt%ntm[j,0])+str(fmt%ntm[j,1])+str(fmt%ntm[j,2]))
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

