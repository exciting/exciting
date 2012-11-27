
import numpy as np

from vtk import vtkPointData, vtkDataArray, vtkUnstructuredGrid, vtkPoints, \
                vtkIdList, vtkStructuredPoints
from ase.visualize.vtk.cell import vtkUnitCellModule
from ase.visualize.vtk.data import vtkDataArrayFromNumPyBuffer, \
                                   vtkDoubleArrayFromNumPyArray, \
                                   vtkDoubleArrayFromNumPyMultiArray

# -------------------------------------------------------------------

class vtkBaseGrid:
    def __init__(self, npoints, cell):
        self.npoints = npoints

        # Make sure cell argument is correct type
        assert isinstance(cell, vtkUnitCellModule)
        self.cell = cell

        self.vtk_pointdata = None

    def set_point_data(self, vtk_pointdata):
        if self.vtk_pointdata is not None:
            raise RuntimeError('VTK point data already present.')

        assert isinstance(vtk_pointdata, vtkPointData)
        self.vtk_pointdata = vtk_pointdata
        #self.vtk_pointdata.SetCopyScalars(False)
        #self.vtk_pointdata.SetCopyVectors(False)
        #self.vtk_pointdata.SetCopyNormals(False)

    def get_point_data(self):
        if self.vtk_pointdata is None:
            raise RuntimeError('VTK point data missing.')

        return self.vtk_pointdata

    def get_number_of_points(self):
        return self.npoints

    def add_scalar_data_array(self, data, name=None, active=True):

        # Are we converting from NumPy buffer to VTK array?
        if isinstance(data, vtkDataArray):
            vtk_sda = data
        elif isinstance(data, vtkDataArrayFromNumPyBuffer):
            vtk_sda = data.get_output()
        else:
            raise ValueError('Data is not a valid scalar data array.')

        del data

        assert vtk_sda.GetNumberOfComponents() == 1
        assert vtk_sda.GetNumberOfTuples() == self.npoints

        if name is not None:
            vtk_sda.SetName(name)

        # Add VTK array to VTK point data
        self.vtk_pointdata.AddArray(vtk_sda)

        if active:
            self.vtk_pointdata.SetActiveScalars(name)

        return vtk_sda

    def add_vector_data_array(self, data, name=None, active=True):

        # Are we converting from NumPy buffer to VTK array?
        if isinstance(data, vtkDataArray):
            vtk_vda = data
        elif isinstance(data, vtkDataArrayFromNumPyBuffer):
            vtk_vda = data.get_output()
        else:
            raise ValueError('Data is not a valid vector data array.')

        del data

        assert vtk_vda.GetNumberOfComponents() == 3
        assert vtk_vda.GetNumberOfTuples() == self.npoints

        if name is not None:
            vtk_vda.SetName(name)

        # Add VTK array to VTK point data
        self.vtk_pointdata.AddArray(vtk_vda)

        if active:
            self.vtk_pointdata.SetActiveVectors(name)

        return vtk_vda

# -------------------------------------------------------------------

class vtkAtomicPositions(vtkBaseGrid):
    """Provides an interface for adding ``Atoms``-centered data to VTK
    modules. Atomic positions, e.g. obtained using atoms.get_positions(),
    constitute an unstructured grid in VTK, to which scalar and vector
    can be added as point data sets.

    Just like ``Atoms``, instances of ``vtkAtomicPositions`` can be divided
    into subsets, which makes it easy to select atoms and add properties.

    Example:

    >>> cell = vtkUnitCellModule(atoms)
    >>> apos = vtkAtomicPositions(atoms.get_positions(), cell)
    >>> apos.add_scalar_property(atoms.get_charges(), 'charges')
    >>> apos.add_vector_property(atoms.get_forces(), 'forces')

    """
    def __init__(self, pos, cell):
        """Construct basic VTK-representation of a set of atomic positions.

        pos: NumPy array of dtype float and shape ``(n,3)``
            Cartesian positions of the atoms.
        cell: Instance of vtkUnitCellModule of subclass thereof
            Holds information equivalent to that of atoms.get_cell().

        """
        # Make sure position argument is a valid array
        if not isinstance(pos, np.ndarray):
            pos = np.array(pos)

        assert pos.dtype == float and pos.shape[1:] == (3,)

        vtkBaseGrid.__init__(self, len(pos), cell)

        # Convert positions to VTK array
        npy2da = vtkDoubleArrayFromNumPyArray(pos)
        vtk_pda = npy2da.get_output()
        del npy2da

        # Transfer atomic positions to VTK points
        self.vtk_pts = vtkPoints()
        self.vtk_pts.SetData(vtk_pda)

        # Create a VTK unstructured grid of these points
        self.vtk_ugd = vtkUnstructuredGrid()
        self.vtk_ugd.SetWholeBoundingBox(self.cell.get_bounding_box())
        self.vtk_ugd.SetPoints(self.vtk_pts)

        # Extract the VTK point data set
        self.set_point_data(self.vtk_ugd.GetPointData())

    def get_points(self, subset=None):
        """Return (subset of) vtkPoints containing atomic positions.

        subset=None: list of int
            A list of indices into the atomic positions; ignored if None.

        """
        if subset is None:
            return self.vtk_pts

        # Create a list of indices from the subset
        vtk_il = vtkIdList()
        for i in subset:
            vtk_il.InsertNextId(i)

        # Allocate VTK points for subset
        vtk_subpts = vtkPoints()
        vtk_subpts.SetDataType(self.vtk_pts.GetDataType())
        vtk_subpts.SetNumberOfPoints(vtk_il.GetNumberOfIds())

        # Transfer subset of VTK points
        self.vtk_pts.GetPoints(vtk_il, vtk_subpts)

        return vtk_subpts

    def get_unstructured_grid(self, subset=None):
        """Return (subset of) an unstructured grid of the atomic positions.

        subset=None: list of int
            A list of indices into the atomic positions; ignored if None.

        """
        if subset is None:
            return self.vtk_ugd

        # Get subset of VTK points
        vtk_subpts = self.get_points(subset)

        # Create a VTK unstructured grid of these points
        vtk_subugd = vtkUnstructuredGrid()
        vtk_subugd.SetWholeBoundingBox(self.cell.get_bounding_box())
        vtk_subugd.SetPoints(vtk_subpts)

        return vtk_subugd

    def add_scalar_property(self, data, name=None, active=True):
        """Add VTK-representation of scalar data at the atomic positions.

        data: NumPy array of dtype float and shape ``(n,)``
            Scalar values corresponding to the atomic positions.
        name=None: str
            Unique identifier for the scalar data.
        active=True: bool
            Flag indicating whether to use as active scalar data.

        """
        # Make sure data argument is a valid array
        if not isinstance(data, np.ndarray):
            data = np.array(data)

        assert data.dtype == float and data.shape == (self.npoints,)

        # Convert scalar properties to VTK array
        npa2da = vtkDoubleArrayFromNumPyArray(data)
        return vtkBaseGrid.add_scalar_data_array(self, npa2da, name, active)

    def add_vector_property(self, data, name=None, active=True):
        """Add VTK-representation of vector data at the atomic positions.

        data: NumPy array of dtype float and shape ``(n,3)``
            Vector components corresponding to the atomic positions.
        name=None: str
            Unique identifier for the vector data.
        active=True: bool
            Flag indicating whether to use as active vector data.

        """
        # Make sure data argument is a valid array
        if not isinstance(data, np.ndarray):
            data = np.array(data)

        assert data.dtype == float and data.shape == (self.npoints,3,)

        # Convert vector properties to VTK array
        npa2da = vtkDoubleArrayFromNumPyArray(data)
        return vtkBaseGrid.add_vector_data_array(self, npa2da, name, active)

# -------------------------------------------------------------------

class vtkVolumeGrid(vtkBaseGrid):
    def __init__(self, elements, cell, origin=None):

        # Make sure element argument is a valid array
        if not isinstance(elements, np.ndarray):
            elements = np.array(elements)

        assert elements.dtype == int and elements.shape == (3,)
        self.elements = elements

        vtkBaseGrid.__init__(self, np.prod(self.elements), cell)

        # Create a VTK grid of structured points
        self.vtk_spts = vtkStructuredPoints()
        self.vtk_spts.SetWholeBoundingBox(self.cell.get_bounding_box())
        self.vtk_spts.SetDimensions(self.elements)
        self.vtk_spts.SetSpacing(self.get_grid_spacing())

        if origin is not None:
            self.vtk_spts.SetOrigin(origin)

        # Extract the VTK point data set
        self.set_point_data(self.vtk_spts.GetPointData())

    def get_grid_spacing(self):
        # Periodic boundary conditions leave out one boundary along an axis
        # Zero/fixed boundary conditions leave out both boundaries of an axis
        return self.cell.get_size()/(self.elements+1.0-self.cell.get_pbc())

    def get_relaxation_factor(self):
        # The relaxation factor is a floating point value between zero and one.
        # It expresses the need for smoothening (relaxation) e.g. of isosurfaces
        # due to coarse grid spacings. Larger grid spacing -> larger relaxation.
        x = self.get_grid_spacing().mean()/self.cell.get_characteristic_length()

        # The relaxation function f(x) satisfies the following requirements
        # f(x) -> 0 for x -> 0+   and   f(x) -> b for x -> inf
        # f'(x) -> a for x -> 0+  and   f'(x) -> 0 for x -> inf

        # Furthermore, it is a rescaling of arctan, hence we know
        # f(x) = 2 b arctan(a pi x / 2 b) / pi

        # Our reference point is x = r for which medium relaxion is needed
        # f(r) = b/2   <=>   r = 2 b / a pi   <=>   a = 2 b / r pi
        r = 0.025 # corresponding to 0.2 Ang grid spacing in 8 Ang cell
        b = 0.5
        f = 2*b*np.arctan(x/r)/np.pi

        if f > 0.1:
           return f.round(1)
        else:
           return None

    def get_structured_points(self):
        return self.vtk_spts

    def add_scalar_field(self, data, name=None, active=True):

        # Make sure data argument is a valid array
        if not isinstance(data, np.ndarray):
            data = np.array(data)

        assert data.dtype == float and data.shape == tuple(self.elements)

        # Convert scalar field to VTK array
        npa2da = vtkDoubleArrayFromNumPyMultiArray(data[...,np.newaxis])
        return vtkBaseGrid.add_scalar_data_array(self, npa2da, name, active)

    def add_vector_field(self, data, name=None, active=True):

        # Make sure data argument is a valid array
        if not isinstance(data, np.ndarray):
            data = np.array(data)

        assert data.dtype == float and data.shape == tuple(self.elements)+(3,)

        # Convert vector field to VTK array
        npa2da = vtkDoubleArrayFromNumPyMultiArray(data)
        return vtkBaseGrid.add_vector_data_array(self, npa2da, name, active)

