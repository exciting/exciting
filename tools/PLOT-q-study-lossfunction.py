#!/usr/bin/python
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from lxml import etree
import os
import sys


#-------------------------------------------------------------------------------
#Check arguments
narg=len(sys.argv)-1
if narg!=3:
    print "\nIncorrect number of arguments\n\n*** Usage: ***\n $>PLOT-lossfkt.py EMIN EMAX OFFSET \nThe spectra will be ploted in the energy range [EMIN,EMAX]. \nOFFSET is a numerical value for the offset between the curves.\n**************\n"
    sys.exit()



emin=float(sys.argv[1])
emax=float(sys.argv[2])
offset=float(sys.argv[3])

if emin>=emax:
    print "\nEMIN should be smaller than EMAX. Exiting\n\n*** Usage: ***\n $>PLOT-lossfkt.py EMIN EMAX OFFSET \nThe spectra will be ploted in the energy range [EMIN,EMAX]. \nOFFSET is a numerical value for the offset between the curves.\n**************\n"
    sys.exit()


#-------------------------------------------------------------------------------
#Parse LOSS function data files
xdata=[]
ydata=[]
labels=[]
legends=[]
directory = os.getcwd()
fnames = [filen for filen in os.listdir(directory) if filen.lower().startswith("loss") and filen.lower().endswith("xml") ]
nfiles=len(fnames)
fnames = sorted(fnames,reverse=True)
for i,fname in enumerate(fnames):
    xdata.append([])
    ydata.append([])
    labels.append({})
    print "Parsing "+fname
    sfname = fname.split("/")[-1].split("_")
    if "FXCRPA" in sfname: legend="RPA "
    if "FXC00" in sfname: legend="RPA "
    if "FXCALDA" in sfname: legend=legend+"ALDA "
    if "FXC05" in sfname: legend=legend+"ALDA "
    if not "NLF" in sfname: legend=legend+"(LFE) "
    if "NLF" in sfname: legend=legend+"(no-LFE) "
    if sfname[-2][0:2]=="OC": legend=legend+"Optical(%s)"%(sfname[-2][2:])
    legends.append(legend)
    tree=etree.parse(fname)
    labels[i]["xlabel"] = tree.xpath('/loss/mapdef/variable1')[0].attrib["name"]
    labels[i]["ylabel"] = tree.xpath('/loss/mapdef/function1')[0].attrib["name"]
    
    for elem in tree.xpath('/loss/map'):
        freq = float(elem.attrib["variable1"]) 
        if emin <= freq <= emax:
            xdata[i].append(float(elem.attrib["variable1"]))
            ydata[i].append(float(elem.attrib["function1"])+i*float(offset))

#-------------------------------------------------------------------------------
#Plot LOSS function/s 
fig=plt.figure(1,figsize=(8,5.5))

params = {'font.size':15,
          'xtick.major.size': 5,
          'ytick.major.size': 5,
          'patch.linewidth': 1.5,
          'axes.linewidth': 2.,
          'axes.formatter.limits': (-4, 6),
          'lines.linewidth': 1.8,
          'lines.markeredgewidth':2.0,
          'lines.markersize':18,
          'legend.fontsize':11,
          'legend.borderaxespad':1,
          'legend.borderpad':0.5,
          'savefig.dpi':80}

plt.rcParams.update(params)

ax=fig.add_subplot(111)

for i in range(nfiles):
    ax.plot(xdata[i],ydata[i],'0.3',label=legends[i])

ax.set_xlim(emin,emax)
ax.set_ylim(0)
ax.set_xlabel(str.capitalize(labels[0]["xlabel"])+" [eV]")
ax.set_ylabel(str.capitalize(labels[0]["ylabel"]))

plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png')
plt.show()
