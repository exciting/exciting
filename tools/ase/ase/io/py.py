from ase.atoms import Atoms


def write_py(fileobj, images, **kwargs):
    if isinstance(fileobj, str):
        fileobj = open(fileobj, 'w')

    fileobj.write('from ase import Atoms\n\n')
    fileobj.write('import numpy as np\n\n')
    
    if not isinstance(images, (list, tuple)):
        images = [images]
    fileobj.write('images = [\n')

    for image in images:
        fileobj.write("    Atoms(symbols='%s',\n"
                      "          pbc=np.%s,\n"
                      "          cell=np.array(\n      %s,\n"
                      "          positions=np.array(\n      %s),\n" % (
            image.get_chemical_symbols(reduce=True),
            repr(image.pbc),
            repr(image.cell)[6:],
            repr(image.positions)[6:]))
        
    fileobj.write(']')
