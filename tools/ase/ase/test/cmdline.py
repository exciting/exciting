import numpy as np
from ase.tasks.main import run
atoms, task = run('H2 -R 0.01 -F 5,2 --atomize')
atoms, task = run('H2 H -s')
data = task.data['H2']
assert abs(data['relaxed energy'] - 1.0705) < 0.0001
assert abs(data['distance'] - 0.7790) < 0.0001
assert abs(data['frequency'] - 0.8628) < 0.0001
assert abs(data['atomic energy'] - data['relaxed energy'] - 5.3495) < 0.0001

atoms, task = run('bulk Cu -F 5,2')
atoms, task = run('bulk Cu -s')
data = task.data['Cu']
assert abs(data['fitted energy'] - -0.0070) < 0.0001
assert abs(data['volume'] - 11.5654) < 0.0001
assert abs(data['B'] - 0.83910) < 0.00001
