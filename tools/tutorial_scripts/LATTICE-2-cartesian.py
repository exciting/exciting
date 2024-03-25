#!/usr/bin/python
#
# @Pasquale Pavone  (2011, April, 4)
#_______________________________________________________________________________

from   lxml  import etree
from   sys  import stdin
from   numpy import *
import subprocess
import os.path
import shutil
import numpy 
import glob
import math
import sys
import os

#-------------------------------------------------------------------------------

def inpsure(inpfile):
    if (not(str(os.path.exists(inpfile)))): 
        sys.exit("\n ERROR: file "+inpfile+" not found!\n")
    inpname = open(inpfile,"r")
    return inpname

def suffix(ip):
    index = str(ip)
    if (ip < 10): index = '0'+str(ip)
    return index

def read_espresso(file,xfactor):
    ifile = inpsure(file)
    while True:
        line = ifile.readline()
        if (len(line) == 0): break
        if (line.strip().split()[0] == "celldm(1)"):
            alat = float(line.strip().split()[2]) 
        if (line.strip().split()[0] == "CELL_PARAMETERS"):
            basevect=[]
            for i in range(3):
                l = input_file.readline().strip().split()
                for j in range(3):
                    basevect.append(float(l[j]))
    axis_matrix = reshape(numpy.array(basevect),(3,3)) 
    ifile.close()
    factor=alat*xfactor
    return axis_matrix, factor
    
def read_exciting(file,xfactor):
    ifile = inpsure(file)
    input_doc  = etree.parse(ifile)
    input_rut  = input_doc.getroot()
    alat = map(float,input_doc.xpath('/input/structure/crystal/@scale'))[0]
    lst_basevect = input_doc.xpath('//basevect/text()')
    xml_basevect = []
    for ind_basevect in lst_basevect:
        xml_basevect.append(map(float,ind_basevect.split()))
    axis_matrix = reshape(numpy.array(xml_basevect),(3,3)) 
    ifile.close()
    factor=alat*xfactor
    return axis_matrix,factor
    
def read_vasp(file,xfactor):
    sys.exit("\n ERROR: VASP version not yet implemented!\n")
    return
 
#-------------------------------------------------------------------------------

file_lat   = "opt-relative-geometry-lattice"
file_car   = "opt-relative-geometry-cartesian"

input_lat  = inpsure(file_lat)
output_car = open(file_car,"w")

#-------------------------------------------------------------------------------

lespresso = os.path.exists("quantum-espresso")
lexciting = os.path.exists("exciting")
lvasp     = os.path.exists("vasp")
#-------------------------------------------------------------------------------

if (lexciting): int_car = open("internal-vs-strain","w")

#-------------------------------------------------------------------------------

#list_inp = sorted(glob.glob('input-*'))
list_inp = sorted(glob.glob('rundir-*'))

#list_inp = []
#list_inp.append("input-01.xml")

#-------------------------------------------------------------------------------

bohr=0.529177
fmt = '%16.10f'

for inpr in list_inp:
    inp = "input-"+inpr[-2:]+".xml"
    if (lespresso): axis_matrix,factor = read_espresso(inp,bohr)
    if (lexciting): axis_matrix,factor = read_exciting(inp,bohr)
    if (lvasp):     axis_matrix,factor = read_vasp(inp,1.)
    line = input_lat.readline().strip().split()
    pos_lat = []
    eta = float(line[0])
    for i in range(3): pos_lat.append(float(line[i+1])*factor)
    pos_lat_vector = reshape(numpy.array(pos_lat),(3,1)) 
    pos_car_vector = dot(transpose(axis_matrix),pos_lat_vector)
    
    pcv = pos_car_vector
        
    print >>output_car, (fmt%eta), 
    print >>output_car, (fmt%pcv[0,0]), (fmt%pcv[1,0]), (fmt%pcv[2,0])
    #print (fmt%eta),(fmt%pcv[0,0]),(fmt%pcv[1,0]),(fmt%pcv[2,0])
    if (lexciting):
        print >>int_car, (fmt%eta), 
        print >>int_car, (fmt%pcv[0,0]), (fmt%pcv[1,0]), (fmt%pcv[2,0])

#-------------------------------------------------------------------------------

output_car.close()
if (lexciting): int_car.close()

#-------------------------------------------------------------------------------




