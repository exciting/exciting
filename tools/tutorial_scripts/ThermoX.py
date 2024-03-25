#!/usr/bin/python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

import sys
import os
from   lxml  import etree
import numpy
import shutil
from   pylab import *
from   scipy import interpolate

#-------------------------------------------------------------------------------

def checkfile(file):
    if (str(os.path.exists(file))=='False'):  
        sys.exit("ERROR: File "+file+" not found!\n")
    return open(file,"r")

#-------------------------------------------------------------------------------

def readfvib(f,t,a):
    return t.xpath('/thermodynamicproperties/'+a+'/map/@'+f)

#-------------------------------------------------------------------------------

narg  = len(sys.argv)-1

if (narg<3): 
    print "\nIncorrect number of arguments. **Usage**:\n\n",
    print "ThermoX.py TMIN TMAX NTSTEPS\n"
    print "Temperatures should be given in Kelvin\n"
    sys.exit()
    
tmin = float(sys.argv[1]) ; tmax = float(sys.argv[2]) ; ntpt = int(sys.argv[3])
tstep = (tmax-tmin)/float(ntpt) ; temp = []
for i in range(ntpt+1): temp.append(tmin+i*tstep)

#-------------------------------------------------------------------------------

print
print '------------------------------------------------------------------------'
print 'Linear thermal-expansion coefficient for cubic systems'
print '------------------------------------------------------------------------'
bzero = input("\nEnter value for the bulk modulus [GPa] >>>> ")
print

#-------------------------------------------------------------------------------

work_directory = 'alpha'
if (os.path.exists(work_directory)): shutil.rmtree(work_directory)
os.mkdir(work_directory)

fdir = []

for idir in range(3):
    ddir = raw_input('\n>>>> Directory '+str(idir+1)+' of '+str(3)+' : ')
    fdir.append(ddir)
    os.chdir(ddir)
    os.system("cp thermo.xml ../alpha/thermo.xml."+str(ddir))     
    os.system("cp input.xml ../alpha/input.xml."+str(ddir))     
    os.chdir("../")

os.chdir(work_directory)

#os.system("ls -la")

iname = [] ; tname = []

for idir in range(3):
    iname.append("input.xml."+str(fdir[idir])) 
    tname.append("thermo.xml."+str(fdir[idir]))

ifile = [] ; tfile = [] ; itree = [] ; ttree = [] ; alat = [] ; vol  = []
tttt  = [] ; ffff  = [] ; ssss  = [] ; cccc  = [] ; vvvv = [] ; fvib = []
dvib  = []

for i in range(len(iname)):
    ifile.append(checkfile(iname[i]))
    tfile.append(checkfile(iname[i]))
    itree.append(etree.parse(iname[i]))
    ttree.append(etree.parse(tname[i]))
    alat.append(map(float,itree[i].xpath('/input/structure/crystal/@scale'))[0])
    vol.append(alat[i]**3/4.)
    tttt.append(readfvib("variable1",ttree[i],'vibrationalfreeenergy')) 
    ffff.append(readfvib("function1",ttree[i],'vibrationalfreeenergy'))  
    ssss.append(readfvib("function1",ttree[i],'vibrationalentropy'))  
    cccc.append(readfvib("function1",ttree[i],'heatcapacity'))  
    vvvv.append(readfvib("function1",ttree[i],'vibrationalenergy'))  
    tttt[i].insert(0,"0.") ; ffff[i].insert(0,ffff[i][0])
    splf = interpolate.splrep(tttt[i],ffff[i],s=0)
    fvib.append(interpolate.splev(temp,splf,der=0))
    dvib.append(interpolate.splev(temp,splf,der=1))

outf   = open('linear-thermal-expansion','w')
GPa2au = 3.398827e-5
bzero  = bzero*GPa2au

for i in range(len(tttt[0])): tttt[0][i]=float(tttt[0][i])
for i in range(len(tttt[1])): tttt[1][i]=float(tttt[1][i])
for i in range(len(tttt[2])): tttt[2][i]=float(tttt[2][i])

tcut= min(max(tttt[0]),max(tttt[1]),max(tttt[2]))

for i in range(len(temp)):
    alpha=(dvib[0][i]-dvib[2][i])/(vol[2]-vol[0])/bzero/3
    if (temp[i] < tcut+1e-10): print>>outf, temp[i], alpha

outf.close()

#-------------------------------------------------------------------------------

style   = " o-"
if (ntpt > 40): style = "-"
xlabel  = " \"Temperature [K]\""
ylabel  = " \"alpha(T)   [1/K]\""
options = style+xlabel+ylabel
os.system("PLOT-plot.py linear-thermal-expansion "+options)

#-------------------------------------------------------------------------------

os.chdir("../")
 
#------------------------------------------------------------------------------- 