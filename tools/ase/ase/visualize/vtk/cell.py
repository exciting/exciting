
import numpy as np

from ase import Atoms

from vtk import vtkOutlineSource, vtkAxesActor, vtkProperty2D, vtkTextProperty
from ase.visualize.vtk.module import vtkModule, vtkPolyDataModule

# -------------------------------------------------------------------

class vtkUnitCellModule(vtkPolyDataModule):
    def __init__(self, atoms):

        assert isinstance(atoms, Atoms)
        self.pbc = atoms.get_pbc()

        cell = atoms.get_cell()

        """
        if not isinstance(cell, np.ndarray):
            cell = np.array(cell)

        if cell.shape == (3,):
            cell = np.diag(cell)

        assert cell.dtype == float and cell.shape == (3, 3)
        """

        self.vtk_outline = vtkOutlineSource()

        if (cell - np.diag(cell.diagonal())).any():
            corners = np.empty((8,3), dtype=float)
            # edges = [map(int,[(i-1)%2==0,i%4>=2,i>=4]) for i in range(8)]
            for c,a in enumerate([(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), \
                                  (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)]):
                corners[c] = np.dot(a, cell)
            self.bbox = np.array(zip(np.min(corners, axis=0), \
                                     np.max(corners, axis=0))).ravel()
            self.vtk_outline.SetCorners(corners.ravel())
            self.vtk_outline.SetBoxTypeToOriented()
        else:
            self.bbox = np.array(zip(np.zeros(3),cell.diagonal())).ravel()
        self.vtk_outline.SetBounds(self.bbox)

        vtkPolyDataModule.__init__(self, self.vtk_outline)

    def get_bounding_box(self):
        return self.bbox

    def get_size(self):
        return self.bbox[1::2]-self.bbox[0::2]

    def get_pbc(self):
        return self.pbc

    def get_characteristic_length(self):
        return np.prod(self.get_size())**(1.0/3.0)

# -------------------------------------------------------------------

class vtkAxesModule(vtkModule):
    def __init__(self, cell):

        assert isinstance(cell, vtkUnitCellModule)
        self.cell = cell

        l0 = self.cell.get_characteristic_length()

        # Create VTK axes actor (not really a VTK actor though)
        self.vtk_ax = vtkAxesActor()
        self.vtk_ax.SetTipTypeToCone()
        self.vtk_ax.SetConeRadius(5e-2*l0)
        self.vtk_ax.SetShaftTypeToCylinder()
        self.vtk_ax.SetCylinderRadius(5e-3*l0)

        # Create VTK two-dimensional property
        p2d = vtkProperty2D()
        p2d.SetDisplayLocationToBackground()

        vtkModule.__init__(self, self.vtk_ax, p2d)

        # Create VTK text property and apply to axes
        vtk_textproperty = vtkTextProperty()
        vtk_textproperty.SetFontSize(14)
        vtk_textproperty.SetBold(True)
        vtk_textproperty.SetItalic(True)
        vtk_textproperty.SetShadow(True)
        vtk_textproperty.SetJustificationToRight()
        vtk_textproperty.SetVerticalJustificationToCentered()

        self.set_text_property(vtk_textproperty)

    def set_actor(self, vtk_act):
        assert isinstance(vtk_act, vtkAxesActor) #fix for non-vtkActor actor
        self.vtk_act = vtk_act

    def set_property(self, vtk_property):
        assert isinstance(vtk_property, vtkProperty2D)
        for vtk_cap in [self.vtk_ax.GetXAxisCaptionActor2D(),
                        self.vtk_ax.GetYAxisCaptionActor2D(),
                        self.vtk_ax.GetZAxisCaptionActor2D()]:
            #vtk_cap.ThreeDimensionalLeaderOn()
            #vtk_cap.LeaderOn()
            vtk_cap.SetProperty(vtk_property)
            vtk_txt = vtk_cap.GetTextActor()
            vtk_txt.SetProperty(vtk_property)

    def set_text_property(self, vtk_textproperty, scaled=False):
        assert isinstance(vtk_textproperty, vtkTextProperty)
        for vtk_cap in [self.vtk_ax.GetXAxisCaptionActor2D(),
                        self.vtk_ax.GetYAxisCaptionActor2D(),
                        self.vtk_ax.GetZAxisCaptionActor2D()]:
            vtk_txt = vtk_cap.GetTextActor()
            vtk_txt.SetScaledText(scaled)
            vtk_txt.SetTextProperty(vtk_textproperty)

