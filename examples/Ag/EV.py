import numpy as nm
import libxml2
from array import *
from libxml2 import xmlAttr
import matplotlib.pyplot as plt
# read data
eVdoc = libxml2.parseFile("./eV.xml")
ctxt = eVdoc.xpathNewContext()
Volume=nm.array(map(float,map(xmlAttr.getContent,ctxt.xpathEval("//@volume"))))
totalEnergy=nm.array(map(float,map(xmlAttr.getContent,ctxt.xpathEval("//@totalEnergy"))))
# make quadratic fit
p=nm.polyfit(Volume,totalEnergy,2)
curve=nm.poly1d(p)
# find root of derivative to get minimum
min=nm.roots(nm.polyder(p))
# x values for plotting polynomial
xa = nm.linspace(Volume[0],Volume[-1],100)
#plot
plt.figure(1)
plt.title('Ag Volume')
plt.ylabel(r'total energy in $[Hartree]$')
plt.xlabel(r'volume in $[Bohr]^3$')
plt.plot(xa,curve(xa),'-')
plt.plot(Volume,totalEnergy,'o')
plt.plot(min,curve(min),'o')
plt.annotate('minimum Volume '+str(min[0]), xy=(min,curve(min)), xycoords='data' ,
              xytext=(min-10,curve(min)+0.002) ,  arrowprops=dict(arrowstyle="->"))
print 'minimum Volume '+str(min[0])

plt.savefig('EV.png')
plt.savefig('EV.pdf')
print 'plot saved as EV.png'
plt.show()
