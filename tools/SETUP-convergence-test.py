#!/usr/bin/python
# -*- coding: utf-8 -*-
#
################################################################################
#
#_______________________________________________________________________________

from   lxml import etree as ET
import os      # Gives access to functions of the operating system
import sys
import shutil
 
#-------------------------------------------------------------------------------

min_ngk  = 2 ; max_ngk  = 2
min_rgkm = 4 ; max_rgkm = 4

narg     = len(sys.argv)-1

if (narg<4):
    print "\nIncorrect number of arguments. **Usage**:\n\n",
    print "SETUP-convergence-test.py ",
    print "NGRIDK_min NGRIDK_max RGKMAX_min RGKMAX_max [DIRECTORYNAME]\n"
    sys.exit()

min_ngk  = int(sys.argv[1]) ; max_ngk  = int(sys.argv[2])+1
min_rgkm = int(sys.argv[3]) ; max_rgkm = int(sys.argv[4])+1

#-------------------------------------------------------------------------------

if (str(os.path.exists('input.xml'))=='False'): 
    sys.exit("\nERROR: Input file input.xml not found!\n")
fileobj=open("input.xml","r") 
               # Definition of the object fileobj as the file input.xml
 
doc = ET.parse(fileobj)  
               # Creates representation of contents (necessary for xml files)
               # This loads the entire XML document 
               # (input.xml) into an ElementTree instance
 
root = doc.getroot() 
               # This is essential: 
               # It stores the reference to the root file you work 
               # with the doc. prefix is a reference to the variable 
               # doc defined above.
               # The "root" of a .xml file is a unique identificator 
               # for such .xml file
 
atoms = doc.xpath('//atom') 
               # Get list of elements named "atom"

groundstate = doc.xpath('//groundstate') 
               # xpath selects concrete fields, like elements or
               # attributes, in a .xml file. Actually xpath is a 
               # rather complex programming language. The 
               # double bar means it refers to all the nodes
               # of the tree that are below groundstate in the 
               # tree that is the .xml file
               
work_directory = 'workdir'
if (narg > 4): work_directory = sys.argv[5]
if (os.path.exists(work_directory)): shutil.rmtree(work_directory)
os.mkdir(work_directory)
os.chdir(work_directory)

for ngk in range(min_ngk, max_ngk, 2):    # sweeping in variables ngk and rgkm

    for rgkm in range(min_rgkm, max_rgkm, 1):
  
        stringk=str(ngk)
        stringr=str(rgkm)
        if (ngk  < 10): stringk='0'+str(ngk)
        if (rgkm < 10): stringr='0'+str(rgkm)
        workdir="ngridk_"+stringk+"-rgkmax_"+stringr
        if not os.path.exists(workdir): os.makedirs(workdir)

        os.system('echo '+str(ngk)+' >> '+str(workdir)+'/ngridk')
        os.system('echo '+str(rgkm)+' >> '+str(workdir)+'/rgkmax')

        groundstate[0].set('ngridk',str(ngk)+" "+str(ngk)+" "+str(ngk)) 
               # Give different values to the variable ngridk, 
               # which belongs to the groundstate field.
               # The symbol + means concatenation of strings
 
        groundstate[0].set('rgkmax',str(rgkm)) 
 
        fileobj2=open(workdir+"/input.xml","w")   #Write changes to input.xml
 
        fileobj2.write(ET.tostring(root,
                                   pretty_print=True,
                                   xml_declaration=True,
                                   encoding='UTF-8'))
               # The line above assigns the modified variable (groundstate) 
               # to the field groundstate
               # of the input.xml file. ET.tostring means element tree 
               # to string in xml
 
        fileobj2.close()
  
