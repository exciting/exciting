"""
The following contains a database of small molecules

Data for the G2/97 database are from
Raghavachari, Redfern, and Pople, J. Chem. Phys. Vol. 106, 1063 (1997).
See http://www.cse.anl.gov/Catalysis_and_Energy_Conversion/Computational_Thermochemistry.shtml for the original files.

All numbers are experimental values, except for coordinates, which are
MP2(full)/6-31G(d) optimized geometries (from http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/G2-97.htm)

Atomic species:
ref: Curtiss et al. JCP 106, 1063 (1997).
'Enthalpy' is the experimental enthalpies of formation at 0K
'thermal correction' is the thermal corrections H(298)-H(0)

Molecular species:
ref: Staroverov et al. JCP 119, 12129 (2003)
'Enthalpy' is the experimental enthalpies of formation at 298K
'ZPE' is the zero-point energies
'thermal correction' is the thermal enthalpy corrections H(298K) - H_exp(0K)
ZPE and thermal corrections are estimated from B3LYP geometries and vibrations.

Experimental ionization potentials are from http://srdata.nist.gov/cccbdb/.

For details about G2-1 and G2-2 sets see doi:10.1063/1.477422.
"""

from ase.data.g2_1 import data as data_g2_1
from ase.data.g2_2 import data as data_g2_2

data = data_g2_1.copy()

data.update(data_g2_2)

from ase.data.g2_1 import atom_names as atom_names_g2_1
from ase.data.g2_1 import molecule_names as molecule_names_g2_1
from ase.data.g2_2 import atom_names as atom_names_g2_2
from ase.data.g2_2 import molecule_names as molecule_names_g2_2

atom_names = []
for a in atom_names_g2_1 + atom_names_g2_2:
    if a not in atom_names:
        atom_names.append(a)
molecule_names = molecule_names_g2_1 + molecule_names_g2_2

from ase.data.g2_2 import get_ionization_energy
from ase.data.g2_2 import get_atomization_energy
