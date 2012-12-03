import os
import ase

tests = [
    'N2Cu_relax.py',
    'Cu_bulk.py',
    'CO_Au111.py',
    'H2.py',
    'C5H12.py',
    'nanoparticle.py', # SLOW
#    'C2_Cu100.py', # Extremely slow
]

for test in tests:
    filename = ase.__path__[0] + '/optimize/test/' + test
    execfile(filename, {})

