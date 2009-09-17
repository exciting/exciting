import numpy as nm
import libxml2
from array import *
from libxml2 import xmlAttr
import matplotlib.pyplot as plt
        
eVdoc = libxml2.parseFile("./eV.xml")
ctxt = eVdoc.xpathNewContext()

Volume=nm.array(map(float,map(xmlAttr.getContent,ctxt.xpathEval("//@volume"))))
totalEnergy=nm.array(map(float,map(xmlAttr.getContent,ctxt.xpathEval("//@totalEnergy"))))

p=nm.polyfit(Volume,totalEnergy,2)
curve=nm.poly1d(p)
xa = nm.linspace(Volume[0],Volume[-1],100)
print xa
print curve(xa)
plt.figure(1)
plt.title('Aluminum Volume')
plt.ylabel(r'total energy in $[Hartree]$')
plt.xlabel(r'volume in $[Bohr]^3$')
plt.plot(xa,curve(xa),'o')
#plt.plot(Volume,totalEnergy,'o')

 
plt.show()