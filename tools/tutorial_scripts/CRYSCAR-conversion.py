#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#  Transform exciting input from cartesian to lattice coordinates and viceversa.
#_______________________________________________________________________________

from   lxml  import etree  as ET
from   numpy import linalg as NL
import numpy               as NU
import sys
import os

#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

if (narg == 0): input_xml = "input.xml"
if (narg == 1): input_xml = str(sys.argv[1])

root = input_xml[0:len(input_xml)-4]

lxmlinput = os.path.exists(input_xml)

if not(lxmlinput): sys.exit("ERROR: exciting input file "+input_xml+" not found!\n")

#-------------------------------------------------------------------------------

input_obj = open(input_xml,"r")
input_doc = ET.parse(input_obj)
input_rut = input_doc.getroot()

inp_scale = 1.0
inp_stretch = [1.,1.,1.]
inp_cartesian = "false"

xml_cartesian = input_doc.xpath('//structure/@cartesian')
if (len(xml_cartesian) > 0): inp_cartesian = xml_cartesian[0]

xml_scale = map(float,input_doc.xpath('//crystal/@scale'))
if (len(xml_scale) > 0): inp_scale = float(xml_scale[0])

xml_stretch = input_doc.xpath('//crystal/@stretch')
if (len(xml_stretch) > 0): 
    xml_stretch[0] = xml_stretch[0].replace('d', 'e')
    xml_stretch[0] = xml_stretch[0].replace('D', 'e')
    inp_stretch = NU.array(map(float,xml_stretch[0].split()))
    
xml_basevect = input_doc.xpath('//basevect/text()')
inp_basevect = []
for i in xml_basevect:
    inp_basevect.append(map(float,i.split()))
axis_matrix = NU.array(inp_basevect) 

#-------------------------------------------------------------------------------

AM = axis_matrix
AM = NU.dot(inp_scale,AM)
for i in range(len(AM)): AM[i] = NU.dot(inp_stretch[i],AM[i])

#-------------------------------------------------------------------------------

fmt = '%16.10f'
amt = '%15.9f'

if (inp_cartesian == "false"):
    for s in input_rut.iter('structure'): s.set('cartesian', 'true')
    outputfl = root+"-car.xml"
    os.system('cp '+input_xml+" "+root+'-lat.xml')
else:
    for s in input_rut.iter('structure'): s.set('cartesian', 'false')
    outputfl = root+"-lat.xml"
    os.system('cp '+input_xml+" "+root+'-car.xml')

xml_atom = input_doc.xpath('//species/atom')
for j in range(len(xml_atom)):   
    atom = NU.array(map(float,xml_atom[j].get("coord").split())) 
    if (inp_cartesian == "true"): atom = NU.dot(NL.inv(NU.transpose(AM)),atom)
    if (inp_cartesian == "false"): atom = NU.dot(NU.transpose(AM),atom)
    line = (amt%atom[0])
    atomline=(str(amt%atom[0])+str(amt%atom[1])+str(amt%atom[2]))
    xml_atom[j].set("coord",atomline)
    
output_obj = open(outputfl,"w")
output_obj.write(ET.tostring(input_rut, method='xml',
                                        pretty_print=True,
                                        xml_declaration=True,
                                        encoding='UTF-8'))
output_obj.close()   

#-------------------------------------------------------------------------------






