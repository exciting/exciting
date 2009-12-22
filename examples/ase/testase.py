import os
from ase import Atoms
from ase.io import *
from ase.calculators import Exciting
from ase.units import Bohr, Hartree
a=Atoms('N3O', [(0, 0, 0), (1, 0, 0), (0, 0, 1),(0.5,0.5 ,0.5)])
write("geo.exi",a)
b=read("geo.exi")

print a
print a.get_positions()
print b
print b.get_positions()

c=Exciting(template = "template.xsl",
           bin = "../../bin/excitingser",
           speciespath = "../../species/"
           )

LiF=read("../LiF/input.xml",format="exi")
LiF.set_calculator(c)
print LiF
print LiF.get_positions()
print LiF.get_scaled_positions()
print "LIF Total Energy:"
print LiF.get_total_energy()



calculator        =Exciting(kpts=(4, 4, 3)  ,
                            maxscl="3",
                            bin = "../../bin/excitingser",
                            speciespath = "../../species"
                            )
templateCalculator=Exciting(template="template.xsl",
                            bin = "../../bin/excitingser",
                            speciespath = "../../species"
                            )
x=3.75*Bohr
Aluminum=Atoms('Al', [(0, 0, 0)] ,cell=[[x,x,0],[x,0.,x],[0.,x,x]])
Aluminum.set_calculator(calculator)
print "Total Energy:"
print Aluminum.get_total_energy()
print "in Hartree"
print Aluminum.get_total_energy() / Hartree
print "Try template"
Aluminum.set_calculator(templateCalculator)
print "Total Energy:"
print Aluminum.get_total_energy()
print "in Hartree"
print Aluminum.get_total_energy() / Hartree

