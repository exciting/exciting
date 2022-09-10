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

lexciting = os.path.exists('exciting')
lespresso = os.path.exists('quantum-espresso')
lvasp     = os.path.exists('vasp')
lsource   = os.path.exists('source.xml')

if   (lexciting):
      if (not(lsource)): sys.exit("\n ERROR: file source.xml not found!\n")
      input_obj = open("source.xml","r")
      input_doc = etree.parse(input_obj)
      input_rut = input_doc.getroot()
      alat      = float(map(float,input_doc.xpath('/input/structure/crystal/@scale'))[0])
      lst_basev = input_doc.xpath('//basevect/text()')
      xml_basev = []
      for ind_basev in lst_basev:
          xml_basev.append(map(float,ind_basev.split()))
      ax_matrix = numpy.array(xml_basev)
      covera    = ax_matrix[2][2]
      
elif (lespresso):
      alat   = 4.6533#input("\nEnter equilibrium lattice constant [bohr] >>>> ")
      covera = 6.0000#input("\nEnter c/a                                 >>>> ")

else: sys.exit("\n ERROR: file source.xml not found!\n")

file_planar = open('planar',"w")
print >> file_planar, alat, covera
file_planar.close()

#-------------------------------------------------------------------------------
