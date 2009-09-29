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
minv=nm.roots(nm.polyder(p))
print 'minimum Volume '+str(minv)
print 'minimum energy at scale '+str(pow(minv/2,1./3.))
# x values for plotting polynomial
xa = nm.linspace(Volume[0],Volume[-1],100)
#plot
plt.figure(1)
plt.title('Ag Volume')
plt.ylabel(r'total energy in $[Hartree]$')
plt.xlabel(r'volume in $[Bohr]^3$')
plt.plot(xa,curve(xa),'-')
plt.plot(Volume,totalEnergy,'o')
plt.annotate('minimum Volume '+str(minv), xy=(minv,curve(minv)), xycoords='data' ,
              xytext=(minv-7,curve(minv)+0.002) ,  arrowprops=dict(arrowstyle="->"))


plt.savefig('EV.png')

print 'plot saved as EV.png'
plt.show()
