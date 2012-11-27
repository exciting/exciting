import os

from ase.test import NotAvailable

try:
    abinit_pp_path = os.getenv('ABINIT_PP_PATH')
    if abinit_pp_path == None:
        raise NotAvailable('ABINIT_PP_PATH not defined')
except NotAvailable:
    raise NotAvailable('Abinit required')

import numpy as np

from ase.tasks.main import run

atoms, task = run("abinit bulk Al -x fcc -a 4.04 --k-point-density=3.0 -p xc='PBE',ecut=340")
atoms, task = run('abinit bulk Al -s')
