
import numpy as np

from vtk import vtkProp3D, vtkPolyDataMapper, vtkActor, vtkLODActor, \
                vtkPointSet, vtkGlyph3D, vtkRenderer
from ase.visualize.vtk.sources import vtkCustomGlyphSource, \
                                      vtkClampedGlyphSource

# -------------------------------------------------------------------

class vtkModule:
    """Modules represent a unified collection of VTK objects needed for
    introducing basic visualization concepts such as surfaces or shapes.

    A common trait of all modules is the need for an actor representation and
    corresponding generic properties such as lighting, color and transparency.

    """ 
    def __init__(self, vtk_act, vtk_property=None):
        """Construct basic VTK-representation of a module.

        vtk_act: Instance of vtkActor or subclass thereof
            A vtkActor represents an entity in a rendering scene.

        vtk_property = None: Instance of vtkProperty or subclass thereof
            A vtkProperty represents e.g. lighting and other surface
            properties of a geometric object, in this case the actor.

        """
        self.vtk_act = None
        self.set_actor(vtk_act)

        if vtk_property is not None:
            self.set_property(vtk_property)

    def set_actor(self, vtk_act):
        """Set the actor representing this module in a rendering scene."""
        assert isinstance(vtk_act, vtkActor)
        self.vtk_act = vtk_act

    def set_property(self, vtk_property):
        """Set the property of the actor representing this module."""
        self.vtk_act.SetProperty(vtk_property)

    def get_actor(self):
        """Return the actor representing this module in a rendering scene."""
        return self.vtk_act

# -------------------------------------------------------------------

class vtkLODModule(vtkModule):
    _vtk_actor_class = vtkLODActor

    def get_lod(self):
        return 100

    def set_actor(self, vtk_act):
        vtkModule.set_actor(self, vtk_act)

        if isinstance(vtk_act, vtkLODActor):
            vtk_act.SetNumberOfCloudPoints(self.get_lod())

class vtkPolyDataModule(vtkModule):
    _vtk_actor_class = vtkActor

    __doc__ = vtkModule.__doc__ + """
    Poly data modules are based on polygonal data sets, which can be mapped
    into graphics primitives suitable for rendering within the VTK framework.

    """
    def __init__(self, vtk_polydata, vtk_property=None):
        """Construct VTK-representation of a module containing polygonals.

        vtk_polydata: Instance of vtkPolyData, subclass thereof or
            vtkPolyDataAlgorithm, which produces vtkPolyData as output.
            A vtkPolyData represents a polygonal data set consisting of
            point and cell attributes, which can be mapped to graphics
            primitives for subsequent rendering.

        vtk_property = None: Instance of vtkProperty or subclass thereof
            A vtkProperty represents e.g. lighting and other surface
            properties of a geometric object, in this case the polydata.

        """
        vtkModule.__init__(self, self._vtk_actor_class(), vtk_property)

        self.vtk_dmap = vtkPolyDataMapper()
        self.vtk_dmap.SetInputConnection(vtk_polydata.GetOutputPort())
        self.vtk_act.SetMapper(self.vtk_dmap)

class vtkGlyphModule(vtkPolyDataModule):
    __doc__ = vtkPolyDataModule.__doc__ + """
    Glyph modules construct these polygonal data sets by replicating a glyph
    source across a specific set of points, using available scalar or vector
    point data to scale and orientate the glyph source if desired.

    Example:

    >>> atoms = molecule('C60')
    >>> cell = vtkUnitCellModule(atoms)
    >>> apos = vtkAtomicPositions(atoms.get_positions(), cell)
    >>> vtk_ugd = apos.get_unstructured_grid()
    >>> glyph_source = vtkAtomSource('C')
    >>> glyph_module = vtkGlyphModule(vtk_ugd, glyph_source)

    """
    def __init__(self, vtk_pointset, glyph_source,
                 scalemode=None, colormode=None):
        """Construct VTK-representation of a module containing glyphs.
        These glyphs share a common source, defining their geometrical
        shape, which is cloned and oriented according to the input data.

        vtk_pointset: Instance of vtkPointSet or subclass thereof
            A vtkPointSet defines a set of positions, which may then be
            assigned specific scalar of vector data across the entire set.

        glyph_source: Instance of ~vtk.vtkCustomGlyphSource or subclass thereof
            Provides the basic shape to distribute over the point set.

        """
        assert isinstance(vtk_pointset, vtkPointSet)
        assert isinstance(glyph_source, vtkCustomGlyphSource)

        # Create VTK Glyph3D based on unstructured grid
        self.vtk_g3d = vtkGlyph3D()
        self.vtk_g3d.SetInput(vtk_pointset)
        self.vtk_g3d.SetSource(glyph_source.get_output())
        self.vtk_g3d.SetScaleFactor(glyph_source.get_scale())

        # Clamping normalizes the glyph scaling to within the range [0,1].
        # Setting a scale factor will then scale these by the given factor.
        if isinstance(glyph_source, vtkClampedGlyphSource):
            self.vtk_g3d.SetClamping(True)
            self.vtk_g3d.SetRange(glyph_source.get_range())

        if scalemode is 'off':
            self.vtk_g3d.SetScaleModeToDataScalingOff()
        elif scalemode is 'scalar':
            self.vtk_g3d.SetScaleModeToScaleByScalar()
        elif scalemode is 'vector':
            self.vtk_g3d.SetScaleModeToScaleByVector()
        elif scalemode is not None:
            raise ValueError('Unrecognized scale mode \'%s\'.' % scalemode)

        if colormode is 'scale':
            self.vtk_g3d.SetColorModeToColorByScale()
        elif colormode is 'scalar':
            self.vtk_g3d.SetColorModeToColorByScalar()
        elif colormode is 'vector':
            self.vtk_g3d.SetColorModeToColorByVector()
        elif colormode is not None:
            raise ValueError('Unrecognized scale mode \'%s\'.' % scalemode)

        vtkPolyDataModule.__init__(self, self.vtk_g3d, glyph_source.get_property())

# -------------------------------------------------------------------

class vtkLabelModule(vtkModule):
    def __init__(self, vtk_pointset, vtk_property=None):

        vtk_Module.__init__(self, vtkActor(), vtk_property)

        assert isinstance(vtk_pointset, vtkPointSet)

        self.vtk_dmap = vtkLabeledDataMapper()
        #self.vtk_dmap.SetLabelModeToLabelIds() #TODO XXX strings!!!
        self.vtk_dmap.GetLabelTextProperty().SetFontSize(12)
        self.vtk_dmap.GetLabelTextProperty().SetJustificationToRight()
        self.vtk_dmap.SetInputConnection(vtk_pointset.GetOutputPort())
        self.vtk_act.SetMapper(self.vtk_dmap)

        #vtk_g3d.GetPointIdsName...

# -------------------------------------------------------------------

class vtkModuleAnchor:
    """
    TODO XXX
    """
    def __init__(self):
        """Construct an anchoring point for various VTK modules.

        """
        self.modules = {}

    def get_module(self, name):
        if name not in self.modules.keys():
            raise RuntimeError('Module \'%s\' does not exists.' % name)

        return self.modules[name]

    def add_module(self, name, module):
        if not isinstance(module, vtkModule):
            raise ValueError('Module must be instance of vtkModule.')

        if name in self.modules.keys():
            raise RuntimeError('Module \'%s\' already exists.' % name)

        self.modules[name] = module

    def get_actor(self, name):
        return self.get_module(name).get_actor()

    def add_actors_to_renderer(self, vtk_renderer, name=None):
        assert isinstance(vtk_renderer, vtkRenderer)
        if name is None:
            for module in self.modules.values():
                vtk_renderer.AddActor(module.get_actor())
        else:
            vtk_renderer.AddActor(self.get_actor(name))

