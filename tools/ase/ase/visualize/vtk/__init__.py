try:
    import vtk
    hasvtk = True
    hasmpi = hasattr(vtk, 'vtkMPIController')
except ImportError:
    hasvtk = False
    hasmpi = False

def requirevtk(code=0, parallel=False):
    from ase.test import NotAvailable
    if not hasvtk:
        # VTK required but not installed, force termination
        # with exit status determined by the code argument.
        raise NotAvailable('VTK is not installed.', code)
    if parallel and not hasmpi:
        # VTK MPI required but not installed, force termination
        # with exit status determined by the code argument.
        raise NotAvailable('VTK is not MPI compatible.', code)

def probe_vtk_kilobyte(default=None):
    if not hasvtk:
        return default

    from vtk import vtkCharArray
    vtk_da = vtkCharArray()
    vtk_da.SetNumberOfComponents(1)
    vtk_da.SetNumberOfTuples(1024**2)

    # Size of 1 MB = 1024**2 bytes in "VTK kilobytes"
    size = vtk_da.GetActualMemorySize()
    if size == 1024:
        return 1024
    elif round(abs(size-1024**2/1e3)) == 0:
        return 1e3
    else:
        return default
