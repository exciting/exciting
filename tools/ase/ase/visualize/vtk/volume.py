
import numpy as np

from vtk import vtkContourFilter, vtkDepthSortPolyData
from ase.visualize.vtk.grid import vtkVolumeGrid
from ase.visualize.vtk.module import vtkPolyDataModule
from ase.visualize.vtk.pipeline import vtkSurfaceSmootherPipeline, \
                                       vtkDepthSortPipeline

# -------------------------------------------------------------------

class vtkIsoSurfaceModule(vtkVolumeGrid, vtkPolyDataModule):
    def __init__(self, data, cell, vtk_renderer, contours=1, depthsort=True):
        #TODO don't require vtk_renderer... implement vtkScene
        #TODO contour values from list or autocontours if int

        # Make sure data argument is a valid array
        if not isinstance(data, np.ndarray):
            data = np.array(data)

        vtkVolumeGrid.__init__(self, data.shape, cell)

        self.vtk_iso = vtkContourFilter() #vtkMarchingContourFilter?
        self.vtk_iso.SetInput(self.get_structured_points()) #TODO non-orthogonal

        self.vtk_iso.SetValue(0, 0.25) 
        self.vtk_iso.SetValue(1, -0.25)

        self.smoothpipe = vtkSurfaceSmootherPipeline(self, vtk_iso)

        #TODO use vtkDepthSortPipeline - but vtkPolyDataModule isn't a pipeline

        self.depthsort = depthsort

        if self.depthsort:
            # The depht sort object is set up to generate scalars representing
            # the sort depth. A camera is assigned for the sorting. The camera
            # defines the sort vector (position and focal point).
            self.vtk_ds = vtkDepthSortPolyData()
            self.vtk_ds.SetCamera(vtk_renderer.GetActiveCamera())
            self.vtk_ds.SetInputConnection(self.vtk_iso.GetOutputPort())
            self.vtk_ds.SetDirectionToBackToFront()
            #vtk_ds.SetVector(1, 1, 1)
            #vtk_ds.SortScalarsOn()
            #vtk_ds.Update()

            vtk_renderer.ResetCamera()

            vtkPolyDataModule.__init__(self, self.vtk_ds)
        else:
            vtkPolyDataModule.__init__(self, self.vtk_iso)

        #TODO add color function
        """
        vtk_dmap = vtk.vtkPolyDataMapper()
        vtk_dmap.ScalarVisibilityOn()
        vtk_dmap.SetScalarRange(0, vtk_sda_x.GetMaxNorm()) #TODO GetMinNorm non-existing!
        vtk_dmap.SetScalarModeToUsePointFieldData()
        vtk_dmap.SelectColorArray("x")
        #vtk_dmap.ColorByArrayComponent("x", 0)
        """

