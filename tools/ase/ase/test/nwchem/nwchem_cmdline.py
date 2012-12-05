import os

from ase.test import NotAvailable

try:
    nwchem_command = os.getenv('NWCHEM_COMMAND')
    if nwchem_command == None:
        raise NotAvailable('NWCHEM_COMMAND not defined')
except NotAvailable:
    raise NotAvailable('Nwchem required')

import numpy as np

from ase.tasks.main import run

atoms, task = run("nwchem molecule O2 O -l -p task=gradient")
atoms, task = run('nwchem molecule O2 O -s')
ae = 2 * task.data['O']['energy'] - task.data['O2']['energy']
assert abs(ae - 6.605) < 1e-3
