import os

from ase.test import NotAvailable

try:
    elk_species_path = os.getenv('ELK_SPECIES_PATH')
    if elk_species_path == None:
        raise NotAvailable('ELK_SPECIES_PATH not defined')
except NotAvailable:
    raise NotAvailable('ELK required')

import numpy as np

from ase.tasks.main import run

atoms, task = run("elk bulk Al -x fcc -a 4.04 --k-point-density=3.0 -p xc='PBE',rgkmax=5.0,tforce=True")
atoms, task = run('elk bulk Al -s')
