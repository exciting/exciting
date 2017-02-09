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

maximum_strain = input("\nEnter maximum Lagrangian strain [smax] >>>> ")
if (maximum_strain == 0): 
    work_directory = 'workdir'
    if (os.path.exists(work_directory)): shutil.rmtree(work_directory)
    os.mkdir(work_directory)
    os.system("cp input.xml workdir/input-01.xml")
    output_info = open('workdir/INFO-elastic-constants',"w")
    print >>output_info,\
        "\nMaximum Lagrangian strain       = 0",\
        "\nNumber of strain values         = 1",\
        "\nVolume of equilibrium unit cell = 1.0 [a.u]^3",\
        "\nDeformation code                = 000000",\
        "\nDeformation label               = single\n"
    output_info.close()
    output_eta = open('workdir/strain-01',"w")
    print >>output_eta, "0.00"
    output_eta.close()
    print
    sys.exit("Single unstrained calculation\n")

if (1 < maximum_strain or maximum_strain < 0): 
    sys.exit("ERROR: Maximum Lagrangian strain is out of range [0-1]!\n")
strain_points = \
              input("\nEnter the number of strain values in [-smax,smax] >>>> ")
strain_points = int(abs(strain_points))
if (3 > strain_points or strain_points > 99): 
    sys.exit("ERROR: Number of strain values is out of range [3-99]!\n")

print
print "------------------------------------------------------------------------"
print " List of deformation codes for strains in Voigt notation"
print "------------------------------------------------------------------------"
print "  0 =>  ( eta,  eta,  eta,    0,    0,    0)  | volume strain "
print "  1 =>  ( eta,    0,    0,    0,    0,    0)  | linear strain along x "
print "  2 =>  (   0,  eta,    0,    0,    0,    0)  | linear strain along y "
print "  3 =>  (   0,    0,  eta,    0,    0,    0)  | linear strain along z "
print "  4 =>  (   0,    0,    0,  eta,    0,    0)  | yz shear strain"
print "  5 =>  (   0,    0,    0,    0,  eta,    0)  | xz shear strain"
print "  6 =>  (   0,    0,    0,    0,    0,  eta)  | xy shear strain"
print "  7 =>  (   0,    0,    0,  eta,  eta,  eta)  | shear strain along (111)"
print "  8 =>  ( eta,  eta,    0,    0,    0,    0)  | xy in-plane strain "
print "  9 =>  ( eta, -eta,    0,    0,    0,    0)  | xy in-plane shear strain"
print " 10 =>  ( eta,  eta,  eta,  eta,  eta,  eta)  | global strain" 
print " 11 =>  ( eta,    0,    0,  eta,    0,    0)  | mixed strain" 
print " 12 =>  ( eta,    0,    0,    0,  eta,    0)  | mixed strain"  
print " 13 =>  ( eta,    0,    0,    0,    0,  eta)  | mixed strain"  
print " 14 =>  ( eta,  eta,    0,  eta,    0,    0)  | mixed strain"  
print "------------------------------------------------------------------------"

deformation_code = input("\nEnter deformation code >>>> ")
if (0 > deformation_code or deformation_code > 14): 
    sys.exit("ERROR: Deformation code is out of range [0-14]!\n")

if (deformation_code == 0 ): dc='EEE000'
if (deformation_code == 1 ): dc='E00000'
if (deformation_code == 2 ): dc='0E0000'
if (deformation_code == 3 ): dc='00E000'
if (deformation_code == 4 ): dc='000E00'
if (deformation_code == 5 ): dc='0000E0'
if (deformation_code == 6 ): dc='00000E'
if (deformation_code == 7 ): dc='000EEE'
if (deformation_code == 8 ): dc='EE0000'
if (deformation_code == 9 ): dc='Ee0000'
if (deformation_code == 10): dc='EEEEEE'
if (deformation_code == 11): dc='E00E00'
if (deformation_code == 12): dc='E000E0'
if (deformation_code == 13): dc='E0000E'
if (deformation_code == 14): dc='EE0E00'

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

output_info = open('INFO-elastic-constants',"w")
print >>output_info, "\nMaximum Lagrangian strain       = ", maximum_strain,\
                     "\nNumber of strain values         = ", strain_points,\
                     "\nVolume of equilibrium unit cell = ", volume, "[a.u]^3",\
                     "\nDeformation code                = ", deformation_code,\
                     "\nDeformation label               = ", dc, "\n"
output_info.close()

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

    if (i+1 < 10): strainfile = 'strain-0'+str(i+1)
    else: strainfile = 'strain-'+str(i+1)
    output_str = open(strainfile,"w")
    fmt = '%11.8f'
    print >>output_str, (fmt%eta)
    output_str.close()

    if (abs(eta) < 0.000001): eta=0.000001
    ep=eta
    if (eta < 0.0): em=abs(eta)
    else: em=-eta

#-------------------------------------------------------------------------------

    e=[]
    for j in range(6):
        ev=0
        if  (dc[j:j+1] == 'E' ): ev=ep
        elif(dc[j:j+1] == 'e' ): ev=em
        elif(dc[j:j+1] == '0' ): ev=0
        else: print "==> ", dc; sys.exit("ERROR: deformation code not allowed!") 
        e.append(ev) 
 
#-------------------------------------------------------------------------------

    eta_matrix=mat([[ e[0]   , e[5]/2., e[4]/2.],
                    [ e[5]/2., e[1]   , e[3]/2.],  
  	            [ e[4]/2., e[3]/2., e[2]   ]])

    one_matrix=mat([[  1.0,  0.0,  0.0],
                    [  0.0,  1.0,  0.0],
	 	    [  0.0,  0.0,  1.0]])

#-------------------------------------------------------------------------------
		   		
    norma=1.0
    inorma=0
    eps_matrix=eta_matrix

    if (linalg.norm(eta_matrix) > 0.7):sys.exit("ERROR: too large deformation!") 

    while ( norma > 1.e-10 ):
        x=eta_matrix-0.5*dot(eps_matrix,eps_matrix)
        norma=linalg.norm(x-eps_matrix)      
        eps_matrix=x
        inorma=inorma+1

    def_matrix=one_matrix+eps_matrix
		
    new_axis_matrix=transpose(dot(def_matrix,transpose(axis_matrix)))
    nam=new_axis_matrix

#-------------------------------------------------------------------------------

    xbv = input_doc.xpath('//crystal/basevect')
    fmt = '%22.16f'
    for j in range(3):
        xbv[j].text = str(fmt%nam[j,0])+str(fmt%nam[j,1])+str(fmt%nam[j,2])+" "
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

