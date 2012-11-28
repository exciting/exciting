from ase.data.g2_1_ref import convert
from ase.data.g2_1_ref import atomization_vasp
from ase.data.g2_1_ref import diatomic

info = {}

info['atomization energy'] = {}

info['atomization energy'].update({'reference': convert(atomization_vasp, 0)})
info['atomization energy'].update({'PBE': convert(atomization_vasp, 1)})
info['atomization energy'].update({'PBE0': convert(atomization_vasp, 3)})

info['bondlength'] = {}

info['bondlength'].update({'reference': convert(diatomic, 0)})
info['bondlength'].update({'PBE': convert(diatomic, 1)})
info['bondlength'].update({'PBE0': convert(diatomic, 3)})
