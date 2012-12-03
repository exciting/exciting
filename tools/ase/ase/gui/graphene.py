# encoding: utf-8
"""nanotube.py - Window for setting up Graphene sheets and ribbons.
"""


import gtk
from gettext import gettext as _
from ase.gui.widgets import pack, cancel_apply_ok, oops
from ase.gui.setupwindow import SetupWindow
from ase.gui.pybutton import PyButton
from ase.gui.status import formula
from ase.structure import graphene_nanoribbon
import ase
import numpy as np

introtext = _("""\
Set up a graphene sheet or a graphene nanoribbon.  A nanoribbon may
optionally be saturated with hydrogen (or another element).\
""")

py_template = """
from ase.structure import nanotube

atoms = nanotube(%(n)i, %(m)i, length=%(length)i, bond=%(bl).3f, symbol=%(symb)s)
"""

label_template = _(""" %(natoms)i atoms: %(symbols)s, Volume: %(volume).3f A<sup>3</sup>""")

class SetupGraphene(SetupWindow):
    "Window for setting up a graphene sheet or nanoribbon."
    def __init__(self, gui):
        SetupWindow.__init__(self)
        self.set_title(_("Graphene"))
        vbox = gtk.VBox()

        # Intoductory text
        self.packtext(vbox, introtext)

        # Choose structure
        label = gtk.Label(_("Structure: "))
        self.struct = gtk.combo_box_new_text()
        for s in (_("Infinite sheet"), _("Unsaturated ribbon"), 
                  _("Saturated ribbon")):
            self.struct.append_text(s)
        self.struct.set_active(0)
    
        pack(vbox, [label, self.struct])

        # Orientation
        label = gtk.Label(_("Orientation: "))
        self.orient = gtk.combo_box_new_text()
        self.orient_text = []
        for s in (_("zigzag"), _("armchair")):
            self.orient.append_text(s)
            self.orient_text.append(s)
        self.orient.set_active(0)
        pack(vbox, [label, self.orient])
        pack(vbox, gtk.Label(""))

        # Choose the element and bond length
        label1 = gtk.Label("Element: ")
        #label.set_alignment(0.0, 0.2)
        self.element = gtk.Entry(max=3)
        self.element.set_text("C")
        self.bondlength = gtk.Adjustment(1.42, 0.0, 1000.0, 0.01)
        label2 = gtk.Label(_("  Bond length: "))
        label3 = gtk.Label(_(u"Å"))
        bond_box = gtk.SpinButton(self.bondlength, 10.0, 3)
        pack(vbox, [label1, self.element, label2, bond_box, label3])

        # Choose the saturation element and bond length
        self.sat_label1 = gtk.Label(_("Saturation: "))
        #label.set_alignment(0.0, 0.2)
        self.element2 = gtk.Entry(max=3)
        self.element2.set_text(_("H"))
        self.bondlength2 = gtk.Adjustment(1.12, 0.0, 1000.0, 0.01)
        self.sat_label2 = gtk.Label(_("  Bond length: "))
        self.sat_label3 = gtk.Label(_(u"Å"))
        self.bond_box = gtk.SpinButton(self.bondlength2, 10.0, 3)
        pack(vbox, [self.sat_label1, self.element2, self.sat_label2,
                    self.bond_box, self.sat_label3])

        self.elementinfo = gtk.Label("")
        self.elementinfo.modify_fg(gtk.STATE_NORMAL,
                                   gtk.gdk.color_parse('#FF0000'))
        pack(vbox, [self.elementinfo])
        pack(vbox, gtk.Label(""))

        # Size
        label1 = gtk.Label(_("Width: "))
        label2 = gtk.Label(_("  Length: "))
        self.n = gtk.Adjustment(1, 1, 100, 1)
        self.m = gtk.Adjustment(1, 1, 100, 1)
        spinn = gtk.SpinButton(self.n, 0, 0)
        spinm = gtk.SpinButton(self.m, 0, 0)
        pack(vbox, [label1, spinn, label2, spinm])

        # Vacuum
        label1 = gtk.Label(_("Vacuum: "))
        self.vacuum = gtk.Adjustment(5.0, 0.0, 1000.0, 0.1)
        label2 = gtk.Label(_(u"Å"))
        vac_box = gtk.SpinButton(self.vacuum, 10.0, 2)
        pack(vbox, [label1, vac_box, label2])
        pack(vbox, gtk.Label(""))

        self.status = gtk.Label("")
        pack(vbox,[self.status])
        pack(vbox,[gtk.Label("")])

        # Buttons
        buts = cancel_apply_ok(cancel=lambda widget: self.destroy(),
                               apply=self.apply,
                               ok=self.ok)
        pack(vbox, [buts], end=True, bottom=True)

        # Finalize setup
        self.makeatoms()
        self.struct.connect('changed', self.makeatoms)
        self.orient.connect('changed', self.makeatoms)
        self.element.connect('activate', self.makeatoms)
        self.bondlength.connect('value-changed', self.makeatoms)
        self.element2.connect('activate', self.makeatoms)
        self.bondlength2.connect('value-changed', self.makeatoms)
        self.n.connect('value-changed', self.makeatoms)
        self.m.connect('value-changed', self.makeatoms)
        self.vacuum.connect('value-changed', self.makeatoms)
        self.update_gui()
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui

    def update_element(self, *args):
        "Called when a new element may have been entered."
        # Assumes the element widget is self.element and that a label
        # for errors is self.elementinfo.  The chemical symbol is
        # placed in self.legalelement - or None if the element is
        # invalid.
        symb = []
        if self.struct.get_active() == 2:
            # Saturated nanoribbon
            elements = (self.element.get_text(), self.element2.get_text())
        else:
            elements = (self.element.get_text(), )
            
        for elem in elements:
            if not elem:
                self.invalid_element(_("  No element specified!"))
                return False
            try:
                z = int(elem)
            except ValueError:
                # Probably a symbol
                try:
                    z = ase.data.atomic_numbers[elem]
                except KeyError:
                    self.invalid_element()
                    return False
            try:
                symb.append(ase.data.chemical_symbols[z])
            except KeyError:
                self.invalid_element()
                return False
        self.elementinfo.set_text("")
        self.legal_element = symb[0]
        if len(symb) == 2:
            self.legal_element2 = symb[1]
        else:
            self.legal_element2 = None
        return True
        
    def update_gui(self, *args):
        # Saturation element is only relevant for saturated nanoribbons
        satur = self.struct.get_active() == 2
        for w in (self.element2, self.bond_box):
            w.set_sensitive(satur)
        # Infinite zigzag sheets must have even width
        if self.struct.get_active() == 0 and self.orient.get_active() == 0:
            if self.n.value % 2 == 1:
                self.n.value += 1
            self.n.lower = 2
            self.n.step_increment = 2
        else:
            self.n.lower = 1
            self.n.step_increment = 1
        
    def makeatoms(self, *args):
        self.update_element()
        self.update_gui()
        if self.legal_element is None or (self.struct.get_active() == 2 and
                                          self.legal_element2 is None):
            self.atoms = None
            self.pybut.python = None
            self.status.set_markup(_("Please specify a consistent set of atoms. "))
        else:
            n = int(self.n.value)
            m = int(self.m.value)
            CC = self.bondlength.value
            vacuum = self.vacuum.value
            orient = self.orient_text[self.orient.get_active()]
            elem = self.legal_element
            if self.struct.get_active() == 0:
                # Extended sheet
                self.atoms = graphene_nanoribbon(n, m, type=orient, C_C=CC,
                                                 vacc=vacuum, sheet=True,
                                                 main_element=elem)
            elif self.struct.get_active() == 1:
                # Unsaturated nanoribbon
                self.atoms = graphene_nanoribbon(n, m, type=orient, C_C=CC,
                                                 vacc=vacuum,
                                                 main_element=elem)
            elif self.struct.get_active() == 2:
                # Saturated nanoribbon
                elem2 = self.legal_element2
                self.atoms = graphene_nanoribbon(n, m, type=orient, C_C=CC,
                                                 C_H=self.bondlength2.value,
                                                 vacuum=vacuum,
                                                 saturated=True,
                                                 main_element=elem,
                                                 saturate_element=elem2)
            else:
                raise RuntimeError("Unknown structure in SetupGraphene!")

        # Now, rotate into the xy plane (ase.gui's default view plane)
        pos = self.atoms.get_positions()
        cell = self.atoms.get_cell()
        pbc = self.atoms.get_pbc()
        cell[1,1], cell[2,2] = cell[2,2], cell[1,1]
        x = pos[:,1].copy()
        z = pos[:,2].copy()
        pos[:,1] = z
        pos[:,2] = x
        self.atoms.set_cell(cell)
        self.atoms.set_positions(pos)
        self.atoms.set_pbc([pbc[0], pbc[2], pbc[1]])
        # Find the heights of the unit cell
        h = np.zeros(3)
        uc = self.atoms.get_cell()
        for i in range(3):
            norm = np.cross(uc[i-1], uc[i-2])
            norm /= np.sqrt(np.dot(norm, norm))
            h[i] = np.abs(np.dot(norm, uc[i]))
        label = label_template % {'natoms'  : self.atoms.get_number_of_atoms(),
                                  'symbols' : formula(self.atoms.get_atomic_numbers()),
                                  'volume'  : self.atoms.get_volume()}
        self.status.set_markup(label)                

    def apply(self, *args):
        self.makeatoms()
        if self.atoms is not None:
            self.gui.new_atoms(self.atoms)
            return True
        else:
            oops(_("No valid atoms."),
                 _("You have not (yet) specified a consistent set of "
                   "parameters."))
            return False

    def ok(self, *args):
        if self.apply():
            self.destroy()
            
            
                                                 
