import numpy as np
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data import atomic_numbers

def get_atoms(cmr_data):
    if type(cmr_data)==str:
        raise RuntimeError('cmr db-file: the specified cmr group file does not contain any images, only references.\n'+
                           'This error could be caused by an older version of CMR - or a group file containing only references to other db-files.')
    positions = cmr_data.get('ase_positions')
    numbers = cmr_data.get('ase_atomic_numbers')
    symbols = cmr_data.get('ase_chemical_symbols')
    cell = cmr_data.get('ase_cell')
    pbc = cmr_data.get('ase_pbc')
    tags = np.array(cmr_data.get('ase_tags'))
    magmoms = np.array(cmr_data.get('ase_magnetic_moments'))
    energy = cmr_data.get('ase_potential_energy')

    forces = cmr_data.get('ase_forces')

    
    if numbers is None and not symbols is None:
        numbers = [atomic_numbers[x] for x in symbols]

    if numbers is None or positions is None:
        raise RuntimeError('cmr db-file: there is no or not enough ase data available in the specified db-file.')
    atoms = Atoms(positions=positions,
                  numbers=numbers,
                  cell=cell,
                  pbc=pbc)

    if tags.any():
        atoms.set_tags(list(tags))

    if magmoms.any():
        atoms.set_initial_magnetic_moments(magmoms)
    else:
        magmoms = None

    atoms.calc = SinglePointCalculator(energy, forces, None, magmoms,
                                           atoms)
    return atoms

def read_db(filename, index):
    import cmr
    r = cmr.read(filename)
    if not r.has_key("ase_positions") and r.is_group():
        hashes = r.get_member_hashes()
        hashes = hashes[index]
        if len(hashes)==0:
            raise RuntimeError('cmr db-file: could not find any group members.\n'+
                               'This error could be caused by an older version of CMR - or a group file containing only references to other db-files.')
        if type(hashes)==list:
            return [get_atoms(r.get_xmldata(hash)) for hash in hashes]
        return get_atoms(r.get_xmldata(hashes))
    else:
        return get_atoms(r)
    
    

def write_db(filename, images, **kwargs):
    try:
        import cmr
        cmr.atoms2cmr(images, **kwargs).write(filename)
    except:
        raise
        raise NotAvailable('CMR version>0.3.2 is required')
    


