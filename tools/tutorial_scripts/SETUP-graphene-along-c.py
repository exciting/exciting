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
    sys.exit(" ERROR: Input file input.xml not found!\n")
    
eta_strain = input("\nEnter lagrangian strain >>>> ")
if not(0.5 > abs(eta_strain)) : 
    sys.exit(" ERROR: Lagrangian strain is out of range [-0.5,0.5]!\n")
eps_strain = -1.+sqrt(1.+2.*eta_strain)
    
maximum_displ = \
           input("\nEnter maximum displacement umax [c/a] >>>> ")
if (1 < maximum_displ or maximum_displ < 0): 
    sys.exit(" ERROR: Maximum displacement is out of range [0-1]!\n")
displ_points = \
              input("\nEnter the number of displacements in [0,umax] >>>> ")
displ_points = int(abs(displ_points))
if (1 > displ_points or displ_points > 99): 
    sys.exit(" ERROR: Number of displacements is out of range [1-99]!\n")

#-------------------------------------------------------------------------------

input_obj = open("input.xml","r")
input_doc = etree.parse(input_obj)
input_rut = input_doc.getroot()
 
xml_scale = map(float,input_doc.xpath('/input/structure/crystal/@scale'))
if (xml_scale == []): 
    sys.exit(" ERROR: There is NO scale attribute in input.xml!\n")
new_scale = float(xml_scale[0]) * (1.+eps_strain)    
    
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

output_info = open('INFO-graphene-along-c',"w")
print >>output_info,\
      "\nLattice parameter (alat)              = ", ('%11.8f'%new_scale), "[a.u.]",\
      "\nMaximum displacement (0,0,umax); umax =  ", maximum_displ, "[c/a]",\
      "\nNumber of displacements               = ", displ_points,\
      "\nLagrangian strain                     = ", '%11.8f'%eta_strain, "\n"
output_info.close()

#-------------------------------------------------------------------------------

delta=displ_points-1
convert=1
if (displ_points <= 1):
    displ_points=1
    convert=-1
    delta=1

displ_step=maximum_displ/delta

#-------------------------------------------------------------------------------

for i in range(0,displ_points):

#-------------------------------------------------------------------------------

    displ=i*displ_step

    if (i+1 < 10): displfile = 'displ-0'+str(i+1)
    else: displfile = 'displ-'+str(i+1)
    output_dsp = open(displfile,"w")
    fmt = '%11.8f'
    print >>output_dsp, (fmt%displ)
    output_dsp.close()

    if (abs(displ) < 0.000001): displ=0.000001
    if (abs(displ-0.50) < 0.000001): displ=0.499999
    
    uuu = mat([[ 0.0, 0.0, 0.0], [ 0, 0, displ]])
    ntm = otm + uuu

#-------------------------------------------------------------------------------

    xsc = input_doc.xpath('//crystal')
    fmt = '%16.10f'
    xsc[0].set("scale",str(fmt%new_scale))
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
                                               xml_declaration=False,
                                               encoding='UTF-8'))
    output_obj.close()

#-------------------------------------------------------------------------------

os.chdir('../')
print 

