import os
from ase import Atoms
from ase.io import read, write
from ase.calculators.exciting import Exciting
from ase.units import Bohr, Hartree
import numpy as np
try:
    import lxml
except ImportError:
    print "not there"
    #raise NotAvailable('This test need lxml module.')

a = Atoms('N3O',
          [(0, 0, 0), (1, 0, 0), (0, 0, 1), (0.5, 0.5, 0.5)],
          pbc=True)
a.new_array("rmt",np.array([-1.,-1.,2.3,2]))
print "rmt repeated cell",a.repeat((1,1,2)).get_array("rmt")
print "rmt ",a.get_array("rmt")
write('geo.exi', a)
b = read('geo.exi')

print "a atoms:",a
print "a pos:",a.get_positions()
print "b atoms:",b
print "b pos:",b.get_positions()

calculator1 = Exciting(dir='excitingtestfiles1',
                      kpts=(4, 4, 3),
                      maxscl=3,
                      title="N3O"
                      #bin='/fshome/chm/git/exciting/bin/excitingser'
                      )
 
calculator3 = Exciting(dir='excitingtestfiles3',
                       paramdict={ "title":{"text()":"N3O"},
                                   "groundstate":{"ngridk":"1 2 3","tforce":"true"},
                                   "structureoptimization":{},
                                   "properties":{"dos":{},
                                   "bandstructure":{"plot1d":{
                                      "path":{ "steps":"100", 
                                        "point":[{"coord":"0.75000   0.50000   0.25000", "label":"W"},    
                                                 {"coord":"0.50000   0.50000   0.50000","label":"L" },   
                                                 {"coord":"0.00000   0.00000   0.00000", "label":"GAMMA"},  
                                                 {"coord":"0.50000   0.50000   0.00000", "label":"X"   },  
                                                 {"coord":"0.75000   0.50000   0.25000", "label":"W"  },  
                                                 {"coord":"0.75000   0.37500   0.37500", "label":"K"  }
                                                 ]}
                                            }
                                        }
                                      }
                                   }
                      )
calculator1.write(a)

calculator3.write(b)