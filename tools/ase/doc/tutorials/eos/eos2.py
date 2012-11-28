from ase.io import read 
from ase.units import kJ
from ase.utils.eos import EquationOfState
configs = read('Ag.traj@0:5')  # read 5 configurations
# Extract volumes and energies:
volumes = [ag.get_volume() for ag in configs]
energies = [ag.get_potential_energy() for ag in configs]
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
print B / kJ * 1.0e24, 'GPa'
eos.plot('Ag-eos.png')
