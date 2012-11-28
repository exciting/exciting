
import numpy as np

from vtk import vtkProperty, vtkSphereSource, vtkArrowSource, vtkConeSource

from ase.data import atomic_numbers
from ase.data import covalent_radii as atomic_radii
from ase.data.colors import jmol_colors as atomic_colors
#from ase.data.colors import cpk_colors as atomic_colors

avg_radius = np.mean(atomic_radii[np.isfinite(atomic_radii)])

# -------------------------------------------------------------------

class vtkCustomGlyphSource:
    def __init__(self, scale, vtk_glyph_source):
        self.scale = scale
        self.vtk_glyph_source = vtk_glyph_source
        self.vtk_property = vtkProperty()

    def get_scale(self):
        return self.scale

    def get_property(self):
        return self.vtk_property

    def get_output(self):
        return self.vtk_glyph_source.GetOutput()

class vtkAtomSource(vtkCustomGlyphSource):
    def __init__(self, name, scale=1, fraction=0.25):
        vtkCustomGlyphSource.__init__(self, scale, vtkSphereSource())

        self.number = atomic_numbers[name]
        self.radius = atomic_radii[self.number]
        self.color = atomic_colors[self.number]

        self.vtk_property.SetColor(self.color[0],self.color[1],self.color[2])
        self.vtk_property.SetInterpolationToPhong()
        self.vtk_property.SetDiffuse(0.7)
        self.vtk_property.SetSpecular(0.4)
        self.vtk_property.SetSpecularPower(20)

        self.vtk_glyph_source.SetPhiResolution(16)
        self.vtk_glyph_source.SetThetaResolution(16)
        self.vtk_glyph_source.SetRadius(fraction*self.radius)

# -------------------------------------------------------------------

class vtkClampedGlyphSource(vtkCustomGlyphSource):
    def __init__(self, scale, vtk_glyph_source, range_min, range_max):
        vtkCustomGlyphSource.__init__(self, scale, vtk_glyph_source)
        self.range = (range_min, range_max,)

    def get_range(self):
        return self.range

class vtkForceSource(vtkClampedGlyphSource):
    def __init__(self, maxnorm, scale=1):
        vtkClampedGlyphSource.__init__(self, scale, vtkArrowSource(),
                                       range_min=0.0, range_max=maxnorm)

        self.vtk_property.SetColor(1.0, 0.25, 0.25) # forces are red
        self.vtk_property.SetInterpolationToPhong()
        self.vtk_property.SetDiffuse(0.7)
        self.vtk_property.SetSpecular(0.4)
        self.vtk_property.SetSpecularPower(20)

        self.vtk_glyph_source.SetShaftResolution(12)
        self.vtk_glyph_source.SetShaftRadius(0.03*avg_radius) #default 0.03
        self.vtk_glyph_source.SetTipResolution(20)
        self.vtk_glyph_source.SetTipLength(0.3*avg_radius) #default 0.35
        self.vtk_glyph_source.SetTipRadius(0.1*avg_radius) #default 0.1

class vtkVelocitySource(vtkClampedGlyphSource):
    def __init__(self, maxnorm, scale=1):
        vtkClampedGlyphSource.__init__(self, scale, vtkConeSource(),
                                       range_min=0.0, range_max=maxnorm)

        self.vtk_property.SetColor(0.25, 0.25, 1.0) # velocities blue
        self.vtk_property.SetInterpolationToPhong()
        self.vtk_property.SetDiffuse(0.9)
        self.vtk_property.SetSpecular(1.0)
        self.vtk_property.SetSpecularPower(50)

        self.vtk_glyph_source.SetAngle(6)
        self.vtk_glyph_source.SetHeight(avg_radius)
        self.vtk_glyph_source.SetResolution(16)
        self.vtk_glyph_source.SetCenter(0.05*avg_radius, 0.0, 0.0)

# -------------------------------------------------------------------

class vtkBondSource(vtkCustomGlyphSource):
    def __init__(self, width, scale=1):
        vtkCustomGlyphSource.__init__(self, scale, vtkCylinderSource())

        self.width = width

        self.vtk_property.SetColor(0.25, 1.0, 0.25) # bonds are green
                                                    # (and so are you)

        self.vtk_glyph_source.SetRadius(self.scale*self.width)
        self.vtk_glyph_source.SetResolution(16)

