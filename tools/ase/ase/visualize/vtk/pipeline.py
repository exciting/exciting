
import numpy as np

from vtk import vtkPolyDataMapper, vtkPolyDataNormals, \
                vtkLinearSubdivisionFilter, vtkSmoothPolyDataFilter, \
                vtkDepthSortPolyData, vtkAlgorithm, vtkAlgorithmOutput
from ase.visualize.vtk.grid import vtkVolumeGrid

# -------------------------------------------------------------------

class vtkPipeline:
    _branchable = False

    def __init__(self, vtkish_data=None):

        self.vtkish_data = None
        self.closed = False
        self.pending = False
        self.filters = tuple()

        if vtkish_data is not None:
            self.set_data(vtkish_data)

    def hasdata(self):
        return (self.vtkish_data is not None)

    def hasfilters(self):
        return len(self.filters)>0

    def set_data(self, vtkish_data):

        assert vtkish_data is not None

        if self.hasdata():
            raise RuntimeError('Pipeline already has input.')
        elif vtkish_data in self:
            raise ValueError('Pipeline loop detected.')

        if isinstance(vtkish_data, vtkPipeline):
            # The embedded pipeline takes over
            vtkish_data.signal_close()

        self.vtkish_data = vtkish_data

        if not self.hasfilters():
            return

        vtkish_inner = self.filters[0]
        
        if isinstance(vtkish_inner, vtkPipeline):
            vtkish_inner.set_data(vtkish_data) #TODO does this work?
        else:
            assert isinstance(vtkish_inner, vtkAlgorithm)

            if isinstance(vtkish_data, vtkPipeline):
                # The embedded pipeline takes over
                vtkish_data.connect(vtkish_inner)
            else:
                assert isinstance(vtkish_data, vtkAlgorithm)
                vtkish_inner.SetInputConnection(vtkish_data.GetOutputPort())

    def isempty(self):
        return not self.hasdata() and not self.hasfilters()

    def isclosed(self):
        if self.pending:
            raise RuntimeError('Pipeline output port state is pending.')
        return self.closed

    def signal_close(self):
        if self.closed:
            raise RuntimeError('Pipeline output port is already closed.')
        elif not self._branchable:
            self.pending = True

    def get_output_port(self):
        if self.closed:
            raise RuntimeError('Pipeline output port is closed.')
        elif self.pending:
            self.closed = True
            self.pending = False

        if self.hasfilters():
            vtkish_outer = self.filters[-1]
        elif self.hasdata():
            vtkish_outer = self.vtkish_data
        else:
            raise RuntimeError('Pipeline output port unavailable.')

        if isinstance(vtkish_outer, vtkPipeline):
            return vtkish_outer.get_output_port()
        else:
            assert isinstance(vtkish_outer, vtkAlgorithm)
            return vtkish_outer.GetOutputPort()

    def set_input_connection(self, vtk_port):
        assert isinstance(vtk_port, vtkAlgorithmOutput)

        vtkish_inner = self.filters[0]

        if isinstance(vtkish_inner, vtkPipeline):
            # Connect must be passed down
            vtkish_inner.set_input_connection(vtk_port)
        else:
            vtkish_inner.SetInputConnection(vtk_port)

    def connect(self, vtkish_filter):
        if vtkish_filter in self:
            raise ValueError('Pipeline loop detected.')

        if isinstance(vtkish_filter, vtkPipeline):
            # Connection must be passed down
            if not self.isempty():
                vtkish_filter.set_input_connection(self.get_output_port())

            # The containing pipeline takes over
            vtkish_filter.signal_close()

        elif not self.isempty():
            assert isinstance(vtkish_filter, vtkAlgorithm)
            vtkish_filter.SetInputConnection(self.get_output_port())

    def append(self, vtkish_filter):
        self.connect(vtkish_filter)
        self.filters += (vtkish_filter,)

    def extend(self, vtkish_filters):
        map(self.append, vtkish_filters)

    def __contains__(self, vtkish_candidate):
        if vtkish_candidate == self:
            return True

        if self.hasdata():
            if isinstance(self.vtkish_data, vtkPipeline):
                if vtkish_candidate in self.vtkish_data:
                    return True
            else:
                if vtkish_candidate == self.vtkish_data:
                    return True

        if vtkish_candidate in self.filters:
            return True

        for vtkish_filter in self.filters:
            if isinstance(vtkish_filter, vtkPipeline) \
                and vtkish_candidate in vtkish_filter:
                return True

        return False

    def __getitem__(self, i): #TODO XXX this is a hack
        return self.filters[i]

# -------------------------------------------------------------------

class vtkPolyDataPipeline(vtkPipeline):
    def __init__(self, vtkish_polydata=None):
        vtkPipeline.__init__(self, vtkish_polydata)

    def connect(self, vtkish_filter):
        if isinstance(vtkish_filter, vtkPolyDataMapper):
            self.signal_close()
        vtkPipeline.connect(self, vtkish_filter)

class vtkSurfaceSmootherPipeline(vtkPolyDataPipeline):
    def __init__(self, grid, vtkish_polydata=None, angle=15):

        vtkPolyDataPipeline.__init__(self, vtkish_polydata)

        # Make sure grid argument is correct type
        assert isinstance(grid, vtkVolumeGrid)
        self.grid = grid

        # Split polys with intersection angles greater than angle
        vtk_dnorm = vtkPolyDataNormals()
        vtk_dnorm.SetFeatureAngle(angle)
        vtk_dnorm.SplittingOn()
        vtk_dnorm.ComputeCellNormalsOff()
        vtk_dnorm.ComputePointNormalsOff()
        self.append(vtk_dnorm)

        relax = self.grid.get_relaxation_factor()

        if relax is not None:
            print 'relax=',relax
            #vtk_subdiv = vtkButterflySubdivisionFilter()
            vtk_subdiv = vtkLinearSubdivisionFilter()
            self.append(vtk_subdiv)
    
            # Smooth out some of the sharp points.
            vtk_smooth = vtkSmoothPolyDataFilter()
            vtk_smooth.SetRelaxationFactor(relax)
            self.append(vtk_smooth)

class vtkDepthSortPipeline(vtkPolyDataPipeline):
    def __init__(self, vtk_renderer, vtkish_polydata=None):

        vtkPolyDataPipeline.__init__(self, vtkish_polydata)

        # The depht sort object is set up to generate scalars representing
        # the sort depth. A camera is assigned for the sorting. The camera
        # defines the sort vector (position and focal point).
        vtk_ds = vtkDepthSortPolyData()
        vtk_ds.SetCamera(vtk_renderer.GetActiveCamera())
        vtk_ds.SetDirectionToBackToFront()
        #vtk_ds.SetVector(1, 1, 1)
        #vtk_ds.SortScalarsOn()
        #vtk_ds.Update()
        self.append(vtk_ds)

        vtk_renderer.ResetCamera()

