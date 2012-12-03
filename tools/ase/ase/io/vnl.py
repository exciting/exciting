import numpy as np

from ase.atoms import Atoms


class VNL:
    def __setstate__(self, data):
        self.data = data

def ac(shape, typecode, data, endian):
    x = np.fromstring(data, typecode)
    try:
        x.shape = shape
    except ValueError:
        x = x[::2].copy()
        x.shape = shape
        
    if np.LittleEndian != endian: 
        return x.byteswap() 
    else: 
        return x 

class VNLUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == 'VNLATKStorage.Core.Sample':
            return VNL
        if name == 'array_constructor':
            return ac
        return pickle.Unpickler.find_class(self, module, name)
    
def read_vnl(filename):
    from cStringIO import StringIO
    vnl = VNLUnpickler(StringIO(ZipFile(filename).read('0_object'))).load()
    conf = vnl.data['__properties__']['Atomic Configuration'].data
    numbers = conf['_dataarray_']
    positions = conf['_positions_'].data['_dataarray_']
    return Atoms(numbers=numbers, positions=positions)
