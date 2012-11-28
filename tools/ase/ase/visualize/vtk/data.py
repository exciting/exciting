
import numpy as np
from numpy.ctypeslib import ctypes

from vtk import vtkDataArray, vtkFloatArray, vtkDoubleArray

if ctypes is None:
    class CTypesEmulator:
        def __init__(self):
            self._SimpleCData = np.number
            self.c_float = np.float32
            self.c_double = np.float64
    try:
        import ctypes
    except ImportError:
        ctypes = CTypesEmulator()

# -------------------------------------------------------------------

class vtkNumPyBuffer:
    def __init__(self, data):
        self.strbuf = data.tostring()
        self.nitems = len(data.flat)

    def __len__(self):
        return self.nitems

    def get_pointer(self):
        # Any C/C++ method that requires a void * can be passed a Python
        # string. No check is done to ensure that the string is the correct
        # size, and the string's reference count is not incremented. Extreme
        # caution should be applied when using this feature.
        return self.strbuf

    def notify(self, obj, event):
        if event == 'DeleteEvent':
            del self.strbuf
        else:
            raise RuntimeError('Event not recognized.')

class vtkDataArrayFromNumPyBuffer:
    def __init__(self, vtk_class, ctype, data=None):

        assert issubclass(ctype, ctypes._SimpleCData)
        self.ctype = ctype

        self.vtk_da = vtk_class()
        assert isinstance(self.vtk_da, vtkDataArray)
        assert self.vtk_da.GetDataTypeSize() == np.nbytes[np.dtype(self.ctype)]

        if data is not None:
            self.read_numpy_array(data)

    def read_numpy_array(self, data):

        if not isinstance(data, np.ndarray):
            data = np.array(data, dtype=self.ctype)

        if data.dtype != self.ctype: # NB: "is not" gets it wrong
            data = data.astype(self.ctype)

        self.vtk_da.SetNumberOfComponents(data.shape[-1])

        # Passing the void* buffer to the C interface does not increase
        # its reference count, hence the buffer is deleted by Python when
        # the reference count of the string from tostring reaches zero.
        # Also, the boolean True tells VTK to save (not delete) the buffer
        # when the VTK data array is deleted - we want Python to do this.
        npybuf = vtkNumPyBuffer(data)
        self.vtk_da.SetVoidArray(npybuf.get_pointer(), len(npybuf), True)
        self.vtk_da.AddObserver('DeleteEvent', npybuf.notify)

    def get_output(self):
        return self.vtk_da

    def copy(self):
        vtk_da_copy = self.vtk_da.NewInstance()
        vtk_da_copy.SetNumberOfComponents(self.vtk_da.GetNumberOfComponents())
        vtk_da_copy.SetNumberOfTuples(self.vtk_da.GetNumberOfTuples())

        assert vtk_da_copy.GetSize() == self.vtk_da.GetSize()

        vtk_da_copy.DeepCopy(self.vtk_da)

        return vtk_da_copy

# -------------------------------------------------------------------

class vtkDataArrayFromNumPyArray(vtkDataArrayFromNumPyBuffer):
    """Class for reading vtkDataArray from 1D or 2D NumPy array.

    This class can be used to generate a vtkDataArray from a NumPy array.
    The NumPy array should be of the form <entries> x <number of components>
    where 'number of components' indicates the number of components in 
    each entry in the vtkDataArray. Note that this form is also expected
    even in the case of only a single component.
    """
    def __init__(self, vtk_class, ctype, data=None, buffered=True):

        self.buffered = buffered

        vtkDataArrayFromNumPyBuffer.__init__(self, vtk_class, ctype, data)

    def read_numpy_array(self, data):
        """Read vtkDataArray from NumPy array"""

        if not isinstance(data, np.ndarray):
            data = np.array(data, dtype=self.ctype)

        if data.dtype != self.ctype: # NB: "is not" gets it wrong
            data = data.astype(self.ctype)

        if data.ndim == 1:
            data = data[:, np.newaxis]
        elif data.ndim != 2:
            raise ValueError('Data must be a 1D or 2D NumPy array.')

        if self.buffered:
            vtkDataArrayFromNumPyBuffer.read_numpy_array(self, data)
        else:
            self.vtk_da.SetNumberOfComponents(data.shape[-1])
            self.vtk_da.SetNumberOfTuples(data.shape[0])

            for i, d_c in enumerate(data):
                for c, d in enumerate(d_c):
                    self.vtk_da.SetComponent(i, c, d)

class vtkFloatArrayFromNumPyArray(vtkDataArrayFromNumPyArray):
    def __init__(self, data):
        vtkDataArrayFromNumPyArray.__init__(self, vtkFloatArray,
                                            ctypes.c_float, data)

class vtkDoubleArrayFromNumPyArray(vtkDataArrayFromNumPyArray):
    def __init__(self, data):
        vtkDataArrayFromNumPyArray.__init__(self, vtkDoubleArray,
                                            ctypes.c_double, data)

# -------------------------------------------------------------------

class vtkDataArrayFromNumPyMultiArray(vtkDataArrayFromNumPyBuffer):
    """Class for reading vtkDataArray from a multi-dimensional NumPy array.

    This class can be used to generate a vtkDataArray from a NumPy array.
    The NumPy array should be of the form <gridsize> x <number of components>
    where 'number of components' indicates the number of components in 
    each gridpoint in the vtkDataArray. Note that this form is also expected
    even in the case of only a single component.
    """
    def __init__(self, vtk_class, ctype, data=None, buffered=True):

        self.buffered = buffered

        vtkDataArrayFromNumPyBuffer.__init__(self, vtk_class, ctype, data)

    def read_numpy_array(self, data):
        """Read vtkDataArray from NumPy array"""

        if not isinstance(data, np.ndarray):
            data = np.array(data, dtype=self.ctype)

        if data.dtype != self.ctype: # NB: "is not" gets it wrong
            data = data.astype(self.ctype)

        if data.ndim <=2:
            raise Warning('This is inefficient for 1D and 2D NumPy arrays. ' +
                          'Use a vtkDataArrayFromNumPyArray subclass instead.')

        if self.buffered:
            # This is less than ideal, but will not copy data (uses views).
            # To get the correct ordering, the grid dimensions have to be
            # transposed without moving the last dimension (the components).
            n = data.ndim-1
            for c in range(n//2):
                data = data.swapaxes(c,n-1-c)

            vtkDataArrayFromNumPyBuffer.read_numpy_array(self, data)
        else:
            self.vtk_da.SetNumberOfComponents(data.shape[-1])
            self.vtk_da.SetNumberOfTuples(np.prod(data.shape[:-1]))

            for c, d_T in enumerate(data.T):
                for i, d in enumerate(d_T.flat):
                    self.vtk_da.SetComponent(i, c, d)

class vtkFloatArrayFromNumPyMultiArray(vtkDataArrayFromNumPyMultiArray):
    def __init__(self, data):
        vtkDataArrayFromNumPyMultiArray.__init__(self, vtkFloatArray,
                                                 ctypes.c_float, data)

class vtkDoubleArrayFromNumPyMultiArray(vtkDataArrayFromNumPyMultiArray):
    def __init__(self, data):
        vtkDataArrayFromNumPyMultiArray.__init__(self, vtkDoubleArray,
                                                 ctypes.c_double, data)

