from ase.atoms import Atoms
from ase.parallel import paropen

def write_findsym(fileobj, images):

    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    symbols = images[0].get_chemical_symbols()
    natoms = len(symbols)

    for atoms in images:
        formula  = atoms.get_chemical_symbols(True)
        accuracy = 1.0e-4

        # Write Comment
        fileobj.write('%s\n' % formula)
        fileobj.write('%f   accuracy\n' % accuracy)
        fileobj.write('1    vectors in cartesian coordinates\n')

        # Write cartesian coordinates of vectors
        for x, y, z in atoms.cell:
            fileobj.write('%22.15f %22.15f %22.15f\n' % (x, y, z))

        fileobj.write('1    no known centering\n')

        fileobj.write('1 0 0 \n')
        fileobj.write('0 1 0 \n')
        fileobj.write('0 0 1 \n')

        fileobj.write('%d\n' % natoms)

        numbers  = atoms.get_atomic_numbers()
        for n in numbers:
            fileobj.write('%d ' % (n))

        fileobj.write('\n')

        for x, y, z in atoms.get_positions():
            fileobj.write('%22.15f %22.15f %22.15f\n' % (x, y, z))

