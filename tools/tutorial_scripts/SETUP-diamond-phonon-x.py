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
print
print "------------------------------------------------------------------------"
print " List of labels for phonon modes at the X point in the diamond structure"
print "------------------------------------------------------------------------"
print "  TA =>  Transverse   acoustic mode"
print "  LA =>  Longitudinal acoustic mode"
print "  TO =>  Transverse   optical  mode"
print "  LO =>  Longitudinal optical  mode"
print "------------------------------------------------------------------------"

phonon_mode = raw_input("\nEnter label for phonon mode [1 choice] >>>> ")
if (phonon_mode != "TA" and\
    phonon_mode != "ta" and\
    phonon_mode != "TO" and\
    phonon_mode != "to" and\
    phonon_mode != "LA" and\
    phonon_mode != "la" and\
    phonon_mode != "LO" and\
    phonon_mode != "lo"): 
    sys.exit("ERROR: Phonon mode not allowed!\n")

#-------------------------------------------------------------------------------

input_obj = open("input.xml","r")
input_doc = etree.parse(input_obj)
input_rut = input_doc.getroot()
 
spe=input_doc.xpath('//species')
atom3=etree.SubElement(spe[0],"atom")
atom4=etree.SubElement(spe[0],"atom")

xml_scale = map(float,input_doc.xpath('/input/structure/crystal/@scale'))
if (xml_scale == []): 
    sys.exit("ERROR: There is NO scale attribute in input.xml!\n")
ref_scale=float(xml_scale[0])

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

xml_ngridk = input_doc.xpath('/input/groundstate/@ngridk')
ref_ngridk = numpy.array(map(int,xml_ngridk[0].split()))

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

otm = mat([[ 0.00, 0.00, 0.00], 
           [ 0.00, 0.50, 0.25],
           [ 0.50, 0.50, 0.50], 
           [ 0.50, 1.00, 0.75]])

#-------------------------------------------------------------------------------

for i in range(0,displ_points):

#-------------------------------------------------------------------------------

    displ=i*displ_step-maximum_displ*convert
    
    if (abs(displ) < 0.000001): displ=0.000001
    
    if (i+1 < 10): displfile = 'displ-0'+str(i+1)
    else: displfile = 'displ-'+str(i+1)
    output_dsp = open(displfile,"w")
    fmt = '%11.8f'
    print >>output_dsp, (fmt%displ)
    output_dsp.close()

    z=1
    if (phonon_mode == "LA" or phonon_mode == "la" ):    
        uuu = mat([[ 0., 0.,  z*displ], 
                   [ 0., 0.,  z*displ],
                   [ 0., 0., -z*displ], 
                   [ 0., 0., -z*displ]])
    if (phonon_mode == "LO" or phonon_mode == "lo" ):    
        uuu = mat([[ 0., 0.,  z*displ], 
                   [ 0., 0., -z*displ],
                   [ 0., 0., -z*displ], 
                   [ 0., 0.,  z*displ]])
    if (phonon_mode == "TA" or phonon_mode == "ta" ):    
        z=sqrt(2.)
        uuu = mat([[ 0.,  z*displ, 0.], 
                   [ 0.,  z*displ, 0.],
                   [ 0., -z*displ, 0.], 
                   [ 0., -z*displ, 0.]])
    if (phonon_mode == "TO" or phonon_mode == "to" ):    
        z=sqrt(2.)
        uuu = mat([[ 0.,  z*displ, 0.], 
                   [ 0., -z*displ, 0.],
                   [ 0., -z*displ, 0.], 
                   [ 0.,  z*displ, 0.]])

    ntm = otm+uuu

#-------------------------------------------------------------------------------

    fmt = '%16.10f'
    ifmt= '%4i'

    xsc = input_doc.xpath('//crystal')
    new_scale=ref_scale/sqrt(2.)
    xsc[0].set("scale",str(fmt%new_scale))

    xbv = input_doc.xpath('//crystal/basevect')
    one=1.
    zero=0.
    root=sqrt(2.)
    xbv[0].text = str(fmt%one)+str(fmt%zero)+str(fmt%zero)+" "
    xbv[1].text = str(fmt%zero)+str(fmt%one)+str(fmt%zero)+" "
    xbv[2].text = str(fmt%zero)+str(fmt%zero)+str(fmt%root)+" "

    xtm=input_doc.xpath('//species/atom')
    for j in range(len(xtm)): 
        xtm[j].\
              set("coord",str(fmt%ntm[j,0])+str(fmt%ntm[j,1])+str(fmt%ntm[j,2]))

    nnk = []
    nkfactor=1./2**(1./6.)
    nnk.append(int(ref_ngridk[0]*nkfactor+0.5))
    nnk.append(int(ref_ngridk[1]*nkfactor+0.5))
    nnk.append(int(ref_ngridk[2]*nkfactor/sqrt(2.)+0.5))
    xnk=input_doc.xpath('//groundstate')
    xnk[0].set("ngridk",str(ifmt%nnk[0])+str(ifmt%nnk[1])+str(ifmt%nnk[2]))

    if (i+1 < 10): outputfile = 'input-0'+str(i+1)+'.xml'
    else: outputfile = 'input-'+str(i+1)+'.xml'

    if (displ> -1e-8):
        output_obj = open(outputfile,"w")
        output_obj.write(etree.tostring(input_rut, method='xml',
                                                   pretty_print=True,
                                                   xml_declaration=False,
                                                   encoding='UTF-8'))
        output_obj.close()    

#-------------------------------------------------------------------------------

output_file = open('X-phonon-calculation',"w")
print >> output_file, "X-phonon-calculation for the", phonon_mode, "mode"
output_file.close()
os.chdir('../')
print 

