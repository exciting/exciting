from ase.data.g2_1_ref import convert
from ase.data.g2_1_ref import atomization

info = {}

info['atomization energy'] = {}

info['atomization energy'].update({'reference': convert(atomization, 0)})
info['atomization energy'].update({'LSD': convert(atomization, 1)})
info['atomization energy'].update({'PBE': convert(atomization, 2)})
info['atomization energy'].update({'RPBE': convert(atomization, 3)})
info['atomization energy'].update({'BLYP': convert(atomization, 4)})
