import os

from ase.test import NotAvailable

try:
    fleur = os.getenv('FLEUR')
    if fleur == None:
        raise NotAvailable('FLEUR not defined')
except NotAvailable:
    raise NotAvailable('Fleur required')

import numpy as np

from ase.tasks.main import run

atoms, task = run("fleur bulk Al -x fcc -a 4.04 --k-point-density=3.0 -p xc=PBE")
atoms, task = run('fleur bulk Al -s')
