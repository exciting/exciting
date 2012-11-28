import numpy as np
#from Numeric import asarray as Numeric_asarray

#from ase.units import Bohr
#from ase.parallel import paropen

fast = False

# -------------------------------------------------------------------

from vtk import vtkStructuredPoints, vtkDoubleArray, vtkXMLImageDataWriter

def write_vti(filename, atoms, data):

    #if isinstance(fileobj, str):
    #    fileobj = paropen(fileobj, 'w')
        
    if isinstance(atoms, list):
        if len(atoms) > 1:
            raise ValueError('Can only write one configuration to a VTI file!')
        atoms = atoms[0]

    if data is None:
        raise ValueError('VTK XML Image Data (VTI) format requires data!')

    data = np.asarray(data)

    if data.dtype == complex:
        data = np.abs(data)

    cell = atoms.get_cell()

    assert np.all(cell==np.diag(cell.diagonal())), 'Unit cell must be orthogonal'

    bbox = np.array(zip(np.zeros(3),cell.diagonal())).ravel()

    # Create a VTK grid of structured points
    spts = vtkStructuredPoints()
    spts.SetWholeBoundingBox(bbox)
    spts.SetDimensions(data.shape)
    spts.SetSpacing(cell.diagonal() / data.shape)
    #spts.SetSpacing(paw.gd.h_c * Bohr)

    #print 'paw.gd.h_c * Bohr=',paw.gd.h_c * Bohr
    #print 'atoms.cell.diagonal() / data.shape=', cell.diagonal()/data.shape
    #assert np.all(paw.gd.h_c * Bohr==cell.diagonal()/data.shape)

    #s = paw.wfs.kpt_u[0].psit_nG[0].copy()
    #data = paw.get_pseudo_wave_function(band=0, kpt=0, spin=0, pad=False)
    #spts.point_data.scalars = data.swapaxes(0,2).flatten()
    #spts.point_data.scalars.name = 'scalars'

    # Allocate a VTK array of type double and copy data
    da = vtkDoubleArray()
    da.SetName('scalars')
    da.SetNumberOfComponents(1)
    da.SetNumberOfTuples(np.prod(data.shape))

    for i,d in enumerate(data.swapaxes(0,2).flatten()):
        da.SetTuple1(i,d)

    # Assign the VTK array as point data of the grid
    spd = spts.GetPointData() # type(spd) is vtkPointData
    spd.SetScalars(da)

    """
    from vtk.util.vtkImageImportFromArray import vtkImageImportFromArray
    iia = vtkImageImportFromArray()
    #iia.SetArray(Numeric_asarray(data.swapaxes(0,2).flatten()))
    iia.SetArray(Numeric_asarray(data))
    ida = iia.GetOutput()
    ipd = ida.GetPointData()
    ipd.SetName('scalars')
    spd.SetScalars(ipd.GetScalars())
    """

    # Save the ImageData dataset to a VTK XML file.
    w = vtkXMLImageDataWriter()

    if fast:
        w.SetDataModeToAppend()
        w.EncodeAppendedDataOff()
    else:
        w.SetDataModeToAscii()

    w.SetFileName(filename)
    w.SetInput(spts)
    w.Write()

# -------------------------------------------------------------------

from vtk import vtkStructuredGrid, vtkPoints, vtkXMLStructuredGridWriter

def write_vts(filename, atoms, data=None):
    raise NotImplementedError

# -------------------------------------------------------------------

from vtk import vtkUnstructuredGrid, vtkPoints, vtkXMLUnstructuredGridWriter

def write_vtu(filename, atoms, data=None):

    #if isinstance(fileobj, str):
    #    fileobj = paropen(fileobj, 'w')
        
    if isinstance(atoms, list):
        if len(atoms) > 1:
            raise ValueError('Can only write one configuration to a VTI file!')
        atoms = atoms[0]

    """
    if data is None:
        raise ValueError('VTK XML Unstructured Grid (VTI) format requires data!')

    data = np.asarray(data)

    if data.dtype == complex:
        data = np.abs(data)
    """

    cell = atoms.get_cell()

    assert np.all(cell==np.diag(cell.diagonal())), 'Unit cell must be orthogonal' #TODO bounding box with general unit cell?!

    bbox = np.array(zip(np.zeros(3),cell.diagonal())).ravel()

    # Create a VTK grid of structured points
    ugd = vtkUnstructuredGrid()
    ugd.SetWholeBoundingBox(bbox)

    """
    # Allocate a VTK array of type double and copy data
    da = vtkDoubleArray()
    da.SetName('scalars')
    da.SetNumberOfComponents(3)
    da.SetNumberOfTuples(len(atoms))

    for i,pos in enumerate(atoms.get_positions()):
        da.SetTuple3(i,pos[0],pos[1],pos[2])
    """
    p = vtkPoints()
    p.SetNumberOfPoints(len(atoms))
    p.SetDataTypeToDouble()
    for i,pos in enumerate(atoms.get_positions()):
        p.InsertPoint(i,pos[0],pos[1],pos[2])


    ugd.SetPoints(p)

    # Assign the VTK array as point data of the grid
    #upd = ugd.GetPointData() # type(spd) is vtkPointData
    #upd.SetScalars(da)


    # Save the UnstructuredGrid dataset to a VTK XML file.
    w = vtkXMLUnstructuredGridWriter()

    if fast:
        w.SetDataModeToAppend()
        w.EncodeAppendedDataOff()
    else:
        w.GetCompressor().SetCompressionLevel(0)
        w.SetDataModeToAscii()

    w.SetFileName(filename)
    w.SetInput(ugd)
    w.Write()


# -------------------------------------------------------------------

def read_vti(filename):
    raise NotImplementedError

def read_vts(filename):
    raise NotImplementedError

def read_vtu(filename):
    raise NotImplementedError

# -------------------------------------------------------------------

from vtk import vtkXMLFileReadTester

def probe_vtkxml(filename):
    """something..."""

    r = vtkXMLFileReadTester()
    r.SetFileName(filename)

    if r.TestReadFile():
        return r.GetFileDataType()
    else:
        return None

