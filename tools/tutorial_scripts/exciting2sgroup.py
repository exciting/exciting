#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   lxml  import etree  as ET
from   numpy import linalg as NL
import numpy               as NU
import sys
import os

#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

if (narg < 2): 
    print "\nIncorrect number of arguments. **Usage**:\n\n",
    print "exciting2sgroup.py excitingINPUTFILE.xml sgroup.in \n"
    sys.exit()

xmlinput = str(sys.argv[1])
sgroupin = str(sys.argv[2])

outputfl = open(sgroupin,"w")

#-------------------------------------------------------------------------------

lxmlinput = os.path.exists(xmlinput)

if not(lxmlinput): sys.exit("ERROR: exciting input file not found!\n")

#-------------------------------------------------------------------------------

input_obj = open(xmlinput,"r")
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

print >>outputfl, "P"

a = NL.norm(AM[0]) ; avec = AM[0]
b = NL.norm(AM[1]) ; bvec = AM[1]
c = NL.norm(AM[2]) ; cvec = AM[2]

alpha = NU.arccos(NU.dot(bvec,cvec)/b/c)/NU.pi*180.
beta  = NU.arccos(NU.dot(avec,cvec)/a/c)/NU.pi*180.
gamma = NU.arccos(NU.dot(avec,bvec)/a/b)/NU.pi*180.

line = (fmt%a)+(fmt%b)+(fmt%c)+(fmt%alpha)+(fmt%beta)+(fmt%gamma)
for i in range(len(line.split())): print >>outputfl, line.split()[i],"",
print >>outputfl

print >>outputfl 
print >>outputfl, len(input_doc.xpath('//species/atom'))

xml_species = input_doc.xpath('//species')

for i in range(len(xml_species)):
    xml_atom = xml_species[i].findall('atom')
    for j in range(len(xml_atom)):   
        atom = NU.array(map(float,xml_atom[j].get("coord").split()))
        if (inp_cartesian == "true"): atom = NU.dot(NL.inv(NU.transpose(AM)),atom)
        line = (amt%atom[0])
        print >>outputfl, line.strip(),
        for k in range(1,len(atom)): print >>outputfl, (amt%atom[k]),
        print >>outputfl
        print >>outputfl, xml_species[i].get("speciesfile")

#-------------------------------------------------------------------------------






