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

fmt = '%22.16f'
amt = '%15.9f'
outputfl = root+"-rot.xml"

#-------------------------------------------------------------------------------

xml_basevect = input_doc.xpath('//basevect/text()')
basevect = []
for i in xml_basevect:
    basevect.append(map(float,i.split()))
am = NU.array(basevect) 

xml_basevect = input_doc.xpath('//basevect')

xml_basevect[0].text = str(fmt%am[1,0])+str(fmt%am[1,1])+str(fmt%am[1,2])+" "
xml_basevect[1].text = str(fmt%am[2,0])+str(fmt%am[2,1])+str(fmt%am[2,2])+" "
xml_basevect[2].text = str(fmt%am[0,0])+str(fmt%am[0,1])+str(fmt%am[0,2])+" "

xml_atom = input_doc.xpath('//species/atom')
for j in range(len(xml_atom)):   
    atom = NU.array(map(float,xml_atom[j].get("coord").split())) 
    atomline=(str(amt%atom[1])+str(amt%atom[2])+str(amt%atom[0]))
    xml_atom[j].set("coord",atomline)
    
#-------------------------------------------------------------------------------
    
output_obj = open(outputfl,"w")
output_obj.write(ET.tostring(input_rut, method='xml',
                                        pretty_print=True,
                                        xml_declaration=False,
                                        encoding='UTF-8'))
output_obj.close()   

#-------------------------------------------------------------------------------






