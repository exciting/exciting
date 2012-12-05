# encoding: utf-8
"""nanoparticle.py - Window for setting up crystalline nanoparticles.
"""

import gtk
from gettext import gettext as _
from copy import copy
from ase.gui.widgets import pack, cancel_apply_ok, oops, help
from ase.gui.setupwindow import SetupWindow
from ase.gui.pybutton import PyButton
import ase
import ase.data
import numpy as np
# Delayed imports:
# ase.cluster.data
from ase.cluster.cubic import FaceCenteredCubic, BodyCenteredCubic, SimpleCubic
from ase.cluster.hexagonal import HexagonalClosedPacked, Graphite
from ase.cluster import wulff_construction

introtext = _("""\
Create a nanoparticle either by specifying the number of layers, or using the
Wulff construction.  Please press the [Help] button for instructions on how to
specify the directions.
WARNING: The Wulff construction currently only works with cubic crystals!
""")

helptext = _("""
The nanoparticle module sets up a nano-particle or a cluster with a given
crystal structure.

1) Select the element, the crystal structure and the lattice constant(s).
   The [Get structure] button will find the data for a given element.

2) Choose if you want to specify the number of layers in each direction, or if
   you want to use the Wulff construction.  In the latter case, you must specify
   surface energies in each direction, and the size of the cluster.

How to specify the directions:
------------------------------

First time a direction appears, it is interpreted as the entire family of
directions, i.e. (0,0,1) also covers (1,0,0), (-1,0,0) etc.  If one of these
directions is specified again, the second specification overrules that specific
direction.  For this reason, the order matters and you can rearrange the
directions with the [Up] and [Down] keys.  You can also add a new direction,
remember to press [Add] or it will not be included.

Example: (1,0,0) (1,1,1), (0,0,1) would specify the {100} family of directions,
the {111} family and then the (001) direction, overruling the value given for
the whole family of directions.
""")

py_template_layers = """
import ase
%(import)s

surfaces = %(surfaces)s
layers = %(layers)s
lc = %(latconst)s
atoms = %(factory)s('%(element)s', surfaces, layers, latticeconstant=lc)

# OPTIONAL: Cast to ase.Atoms object, discarding extra information:
# atoms = ase.Atoms(atoms)
"""

py_template_wulff = """
import ase
from ase.cluster import wulff_construction

surfaces = %(surfaces)s
esurf = %(energies)s
lc = %(latconst)s
size = %(natoms)s  # Number of atoms
atoms = wulff_construction('%(element)s', surfaces, esurf, size, '%(structure)s',
                           rounding='%(rounding)s', latticeconstant=lc)

# OPTIONAL: Cast to ase.Atoms object, discarding extra information:
# atoms = ase.Atoms(atoms)
"""

class SetupNanoparticle(SetupWindow):
    "Window for setting up a nanoparticle."
    # Structures:  Abbreviation, name, 4-index (boolean), two lattice const (bool), factory
    structure_data = (('fcc', _('Face centered cubic (fcc)'), False, False, FaceCenteredCubic),
                      ('bcc', _('Body centered cubic (bcc)'), False, False, BodyCenteredCubic),
                      ('sc',  _('Simple cubic (sc)'), False, False, SimpleCubic),
                      ('hcp', _('Hexagonal closed-packed (hcp)'), True, True, HexagonalClosedPacked),
                      ('graphite', _('Graphite'), True, True, Graphite),
                      )
    #NB:  HCP is broken!
    
    # A list of import statements for the Python window.
    import_names = {'fcc': 'from ase.cluster.cubic import FaceCenteredCubic',
                    'bcc': 'from ase.cluster.cubic import BodyCenteredCubic',
                    'sc': 'from ase.cluster.cubic import SimpleCubic',
                    'hcp': 'from ase.cluster.hexagonal import HexagonalClosedPacked',
                    'graphite': 'from ase.cluster.hexagonal import Graphite',
                    }
    # Default layer specifications for the different structures.
    default_layers = {'fcc': [( (1,0,0), 6),
                              ( (1,1,0), 9),
                              ( (1,1,1), 5)],
                      'bcc': [( (1,0,0), 6),
                              ( (1,1,0), 9),
                              ( (1,1,1), 5)],
                      'sc':  [( (1,0,0), 6),
                              ( (1,1,0), 9),
                              ( (1,1,1), 5)],
                      'hcp': [( (0,0,0,1), 5),
                              ( (1,0,-1,0), 5)],
                      'graphite': [( (0,0,0,1), 5),
                                   ( (1,0,-1,0), 5)]
                      }
    
    def __init__(self, gui):
        SetupWindow.__init__(self)
        self.set_title(_("Nanoparticle"))
        self.atoms = None
        self.no_update = True
        
        vbox = gtk.VBox()

        # Intoductory text
        self.packtext(vbox, introtext)
           
        # Choose the element
        label = gtk.Label(_("Element: "))
        label.set_alignment(0.0, 0.2)
        element = gtk.Entry(max=3)
        self.element = element
        lattice_button = gtk.Button(_("Get structure"))
        lattice_button.connect('clicked', self.set_structure_data)
        self.elementinfo = gtk.Label(" ")
        pack(vbox, [label, element, self.elementinfo, lattice_button], end=True)
        self.element.connect('activate', self.update)
        self.legal_element = False

        # The structure and lattice constant
        label = gtk.Label(_("Structure: "))
        self.structure = gtk.combo_box_new_text()
        self.list_of_structures = []
        self.needs_4index = {}
        self.needs_2lat = {}
        self.factory = {}
        for abbrev, name, n4, c, factory in self.structure_data:
            self.structure.append_text(name)
            self.list_of_structures.append(abbrev)
            self.needs_4index[abbrev] = n4
            self.needs_2lat[abbrev] = c
            self.factory[abbrev] = factory
        self.structure.set_active(0)
        self.fourindex = self.needs_4index[self.list_of_structures[0]]
        self.structure.connect('changed', self.update_structure)
        
        label2 = gtk.Label(_("Lattice constant:  a ="))
        self.lattice_const_a = gtk.Adjustment(3.0, 0.0, 1000.0, 0.01)
        self.lattice_const_c = gtk.Adjustment(5.0, 0.0, 1000.0, 0.01)
        self.lattice_box_a = gtk.SpinButton(self.lattice_const_a, 10.0, 3)
        self.lattice_box_c = gtk.SpinButton(self.lattice_const_c, 10.0, 3)
        self.lattice_box_a.numeric = True
        self.lattice_box_c.numeric = True
        self.lattice_label_c = gtk.Label(" c =")
        pack(vbox, [label, self.structure])
        pack(vbox, [label2, self.lattice_box_a,
                    self.lattice_label_c, self.lattice_box_c])
        self.lattice_label_c.hide()
        self.lattice_box_c.hide()
        self.lattice_const_a.connect('value-changed', self.update)
        self.lattice_const_c.connect('value-changed', self.update)

        # Choose specification method
        label = gtk.Label(_("Method: "))
        self.method = gtk.combo_box_new_text()
        for meth in (_("Layer specification"), _("Wulff construction")):
            self.method.append_text(meth)
        self.method.set_active(0)
        self.method.connect('changed', self.update_gui_method)
        pack(vbox, [label, self.method])
        pack(vbox, gtk.Label(""))
        self.old_structure = None

        frame = gtk.Frame()
        pack(vbox, frame)
        framebox = gtk.VBox()
        frame.add(framebox)
        framebox.show()
        self.layerlabel = gtk.Label("Missing text")  # Filled in later
        pack(framebox, [self.layerlabel])
        # This box will contain a single table that is replaced when
        # the list of directions is changed.
        self.direction_table_box = gtk.VBox()
        pack(framebox, self.direction_table_box)
        pack(self.direction_table_box, 
             gtk.Label(_("Dummy placeholder object")))
        pack(framebox, gtk.Label(""))
        pack(framebox, [gtk.Label(_("Add new direction:"))])
        self.newdir_label = []
        self.newdir_box = []
        self.newdir_index = []
        packlist = []
        for txt in ('(', ', ', ', ', ', '):
            self.newdir_label.append(gtk.Label(txt))
            adj = gtk.Adjustment(0, -100, 100, 1)
            self.newdir_box.append(gtk.SpinButton(adj, 1, 0))
            self.newdir_index.append(adj)
            packlist.append(self.newdir_label[-1])
            packlist.append(self.newdir_box[-1])
        self.newdir_layers = gtk.Adjustment(5, 0, 100, 1)
        self.newdir_layers_box = gtk.SpinButton(self.newdir_layers, 1, 0)
        self.newdir_esurf = gtk.Adjustment(1.0, 0, 1000.0, 0.1)
        self.newdir_esurf_box = gtk.SpinButton(self.newdir_esurf, 10, 3)
        addbutton = gtk.Button(_("Add"))
        addbutton.connect('clicked', self.row_add)
        packlist.extend([gtk.Label("): "),
                         self.newdir_layers_box,
                         self.newdir_esurf_box,
                         gtk.Label("  "),
                         addbutton])
        pack(framebox, packlist)
        self.defaultbutton = gtk.Button(_("Set all directions to default "
                                          "values"))
        self.defaultbutton.connect('clicked', self.default_direction_table)
        self.default_direction_table()

        # Extra widgets for the Wulff construction
        self.wulffbox = gtk.VBox()
        pack(vbox, self.wulffbox)
        label = gtk.Label(_("Particle size: "))
        self.size_n_radio = gtk.RadioButton(None, _("Number of atoms: "))
        self.size_n_radio.set_active(True)
        self.size_n_adj = gtk.Adjustment(100, 1, 100000, 1)
        self.size_n_spin = gtk.SpinButton(self.size_n_adj, 0, 0)
        self.size_dia_radio = gtk.RadioButton(self.size_n_radio,
                                              _("Volume: "))
        self.size_dia_adj = gtk.Adjustment(1.0, 0, 100.0, 0.1)
        self.size_dia_spin = gtk.SpinButton(self.size_dia_adj, 10.0, 2)
        pack(self.wulffbox, [label, self.size_n_radio, self.size_n_spin,
                    gtk.Label("   "), self.size_dia_radio, self.size_dia_spin,
                    gtk.Label(_(u"Å³"))])
        self.size_n_radio.connect("toggled", self.update_gui_size)
        self.size_dia_radio.connect("toggled", self.update_gui_size)
        self.size_n_adj.connect("value-changed", self.update_size_n)
        self.size_dia_adj.connect("value-changed", self.update_size_dia)
        label = gtk.Label(_("Rounding: If exact size is not possible, "
                            "choose the size"))
        pack(self.wulffbox, [label])
        self.round_above = gtk.RadioButton(None, _("above  "))
        self.round_below = gtk.RadioButton(self.round_above, _("below  "))
        self.round_closest = gtk.RadioButton(self.round_above, _("closest  "))
        self.round_closest.set_active(True)
        butbox = gtk.HButtonBox()
        self.smaller_button = gtk.Button(_("Smaller"))
        self.larger_button = gtk.Button(_("Larger"))
        self.smaller_button.connect('clicked', self.wulff_smaller)
        self.larger_button.connect('clicked', self.wulff_larger)
        pack(butbox, [self.smaller_button, self.larger_button])
        buts = [self.round_above, self.round_below, self.round_closest]
        for b in buts:
            b.connect("toggled", self.update)
        buts.append(butbox)
        pack(self.wulffbox, buts, end=True)

        # Information
        pack(vbox, gtk.Label(""))
        infobox = gtk.VBox()
        label1 = gtk.Label(_("Number of atoms: "))
        self.natoms_label = gtk.Label("-")
        label2 = gtk.Label(_("   Approx. diameter: "))
        self.dia1_label = gtk.Label("-")
        pack(infobox, [label1, self.natoms_label, label2, self.dia1_label])
        pack(infobox, gtk.Label(""))
        infoframe = gtk.Frame(_("Information about the created cluster:"))
        infoframe.add(infobox)
        infobox.show()
        pack(vbox, infoframe)
        
        # Buttons
        self.pybut = PyButton(_("Creating a nanoparticle."))
        self.pybut.connect('clicked', self.makeatoms)
        helpbut = help(helptext)
        buts = cancel_apply_ok(cancel=lambda widget: self.destroy(),
                               apply=self.apply,
                               ok=self.ok)
        pack(vbox, [self.pybut, helpbut, buts], end=True, bottom=True)
        self.auto = gtk.CheckButton(_("Automatic Apply"))
        fr = gtk.Frame()
        fr.add(self.auto)
        fr.show_all()
        pack(vbox, [fr], end=True, bottom=True)
        
        # Finalize setup
        self.update_structure()
        self.update_gui_method()
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui
        self.no_update = False

    def default_direction_table(self, widget=None):
        "Set default directions and values for the current crystal structure."
        self.direction_table = []
        struct = self.get_structure()
        for direction, layers in self.default_layers[struct]:
            adj1 = gtk.Adjustment(layers, -100, 100, 1)
            adj2 = gtk.Adjustment(1.0, -1000.0, 1000.0, 0.1)
            adj1.connect("value-changed", self.update)
            adj2.connect("value-changed", self.update)
            self.direction_table.append([direction, adj1, adj2])
        self.update_direction_table()

    def update_direction_table(self):
        "Update the part of the GUI containing the table of directions."
        #Discard old table
        oldwidgets = self.direction_table_box.get_children()
        assert len(oldwidgets) == 1
        oldwidgets[0].hide()
        self.direction_table_box.remove(oldwidgets[0])
        del oldwidgets  # It should now be gone
        tbl = gtk.Table(len(self.direction_table)+1, 7)
        pack(self.direction_table_box, [tbl])
        for i, data in enumerate(self.direction_table):
            tbl.attach(gtk.Label("%s: " % (str(data[0]),)),
                       0, 1, i, i+1)
            if self.method.get_active():
                # Wulff construction
                spin = gtk.SpinButton(data[2], 1.0, 3)
            else:
                # Layers
                spin = gtk.SpinButton(data[1], 1, 0)
            tbl.attach(spin, 1, 2, i, i+1)
            tbl.attach(gtk.Label("   "), 2, 3, i, i+1)
            but = gtk.Button(_("Up"))
            but.connect("clicked", self.row_swap_next, i-1)
            if i == 0:
                but.set_sensitive(False)
            tbl.attach(but, 3, 4, i, i+1)
            but = gtk.Button(_("Down"))
            but.connect("clicked", self.row_swap_next, i)
            if i == len(self.direction_table)-1:
                but.set_sensitive(False)
            tbl.attach(but, 4, 5, i, i+1)
            but = gtk.Button(_("Delete"))
            but.connect("clicked", self.row_delete, i)
            if len(self.direction_table) == 1:
                but.set_sensitive(False)
            tbl.attach(but, 5, 6, i, i+1)
        tbl.show_all()
        self.update()

    def get_structure(self):
        "Returns the crystal structure chosen by the user."
        return self.list_of_structures[self.structure.get_active()]

    def update_structure(self, widget=None):
        "Called when the user changes the structure."
        s = self.get_structure()
        if s != self.old_structure:
            old4 = self.fourindex
            self.fourindex = self.needs_4index[s]
            if self.fourindex != old4:
                # The table of directions is invalid.
                self.default_direction_table()
            self.old_structure = s
            if self.needs_2lat[s]:
                self.lattice_label_c.show()
                self.lattice_box_c.show()
            else:
                self.lattice_label_c.hide()
                self.lattice_box_c.hide()
            if self.fourindex:
                self.newdir_label[3].show()
                self.newdir_box[3].show()
            else:
                self.newdir_label[3].hide()
                self.newdir_box[3].hide()
        self.update()

    def update_gui_method(self, widget=None):
        "Switch between layer specification and Wulff construction."
        self.update_direction_table()
        if self.method.get_active():
            self.wulffbox.show()
            self.layerlabel.set_text(_("Surface energies (as energy/area, "
                                       "NOT per atom):"))
            self.newdir_layers_box.hide()
            self.newdir_esurf_box.show()
        else:
            self.wulffbox.hide()
            self.layerlabel.set_text(_("Number of layers:"))
            self.newdir_layers_box.show()
            self.newdir_esurf_box.hide()
        self.update()

    def wulff_smaller(self, widget=None):
        "Make a smaller Wulff construction."
        n = len(self.atoms)
        self.size_n_radio.set_active(True)
        self.size_n_adj.value = n-1
        self.round_below.set_active(True)
        self.apply()

    def wulff_larger(self, widget=None):
        "Make a larger Wulff construction."
        n = len(self.atoms)
        self.size_n_radio.set_active(True)
        self.size_n_adj.value = n+1
        self.round_above.set_active(True)
        self.apply()
    
    def row_add(self, widget=None):
        "Add a row to the list of directions."
        if self.fourindex:
            n = 4
        else:
            n = 3
        idx = tuple( [int(a.value) for a in self.newdir_index[:n]] )
        if not np.array(idx).any():
            oops(_("At least one index must be non-zero"))
            return
        if n == 4 and np.array(idx)[:3].sum() != 0:
            oops(_("Invalid hexagonal indices",
                 "The sum of the first three numbers must be zero"))
            return
        adj1 = gtk.Adjustment(self.newdir_layers.value, -100, 100, 1)
        adj2 = gtk.Adjustment(self.newdir_esurf.value, -1000.0, 1000.0, 0.1)
        adj1.connect("value-changed", self.update)
        adj2.connect("value-changed", self.update)
        self.direction_table.append([idx, adj1, adj2])
        self.update_direction_table()

    def row_delete(self, widget, row):
        del self.direction_table[row]
        self.update_direction_table()

    def row_swap_next(self, widget, row):
        dt = self.direction_table
        dt[row], dt[row+1] = dt[row+1], dt[row]
        self.update_direction_table()
        
    def update_gui_size(self, widget=None):
        "Update gui when the cluster size specification changes."
        self.size_n_spin.set_sensitive(self.size_n_radio.get_active())
        self.size_dia_spin.set_sensitive(self.size_dia_radio.get_active())

    def update_size_n(self, widget=None):
        if not self.size_n_radio.get_active():
            return
        at_vol = self.get_atomic_volume()
        dia = 2.0 * (3 * self.size_n_adj.value * at_vol / (4 * np.pi))**(1.0/3)
        self.size_dia_adj.value = dia
        self.update()

    def update_size_dia(self, widget=None):
        if not self.size_dia_radio.get_active():
            return
        at_vol = self.get_atomic_volume()
        n = round(np.pi / 6 * self.size_dia_adj.value**3 / at_vol)
        self.size_n_adj.value = n
        self.update()
                
    def update(self, *args):
        if self.no_update:
            return
        self.update_element()
        if self.auto.get_active():
            self.makeatoms()
            if self.atoms is not None:
                self.gui.new_atoms(self.atoms)
        else:
            self.clearatoms()
        self.makeinfo()

    def set_structure_data(self, *args):
        "Called when the user presses [Get structure]."
        if not self.update_element():
            oops(_("Invalid element."))
            return
        z = ase.data.atomic_numbers[self.legal_element]
        ref = ase.data.reference_states[z]
        if ref is None:
            structure = None
        else:
            structure = ref['symmetry']
                
        if ref is None or not structure in self.list_of_structures:
            oops(_("Unsupported or unknown structure",
                   "Element = %s,  structure = %s" % (self.legal_element,
                                                      structure)))
            return
        for i, s in enumerate(self.list_of_structures):
            if structure == s:
                self.structure.set_active(i)
        a = ref['a']
        self.lattice_const_a.set_value(a)
        self.fourindex = self.needs_4index[structure]
        if self.fourindex:
            try:
                c = ref['c']
            except KeyError:
                c = ref['c/a'] * a
            self.lattice_const_c.set_value(c)
            self.lattice_label_c.show()
            self.lattice_box_c.show()
        else:
            self.lattice_label_c.hide()
            self.lattice_box_c.hide()

    def makeatoms(self, *args):
        "Make the atoms according to the current specification."
        if not self.update_element():
            self.clearatoms()
            self.makeinfo()
            return False
        assert self.legal_element is not None
        struct = self.list_of_structures[self.structure.get_active()]
        if self.needs_2lat[struct]:
            # a and c lattice constants
            lc = {'a': self.lattice_const_a.value,
                  'c': self.lattice_const_c.value}
            lc_str = str(lc)
        else:
            lc = self.lattice_const_a.value
            lc_str = "%.5f" % (lc,)
        if self.method.get_active() == 0:
            # Layer-by-layer specification
            surfaces = [x[0] for x in self.direction_table]
            layers = [int(x[1].value) for x in self.direction_table]
            self.atoms = self.factory[struct](self.legal_element, copy(surfaces),
                                              layers, latticeconstant=lc)
            imp = self.import_names[struct]
            self.pybut.python = py_template_layers % {'import': imp,
                                                      'element': self.legal_element,
                                                      'surfaces': str(surfaces),
                                                      'layers': str(layers),
                                                      'latconst': lc_str,
                                                      'factory': imp.split()[-1]
                                                      }
        else:
            # Wulff construction
            assert self.method.get_active() == 1
            surfaces = [x[0] for x in self.direction_table]
            surfaceenergies = [x[2].value for x in self.direction_table]            
            self.update_size_dia()
            if self.round_above.get_active():
                rounding = "above"
            elif self.round_below.get_active():
                rounding = "below"
            elif self.round_closest.get_active():
                rounding = "closest"
            else:
                raise RuntimeError("No rounding!")
            self.atoms = wulff_construction(self.legal_element, surfaces,
                                            surfaceenergies,
                                            self.size_n_adj.value,
                                            self.factory[struct],
                                            rounding, lc)
            self.pybut.python = py_template_wulff % {'element': self.legal_element,
                                                     'surfaces': str(surfaces),
                                                     'energies': str(surfaceenergies),
                                                     'latconst': lc_str,
                                                     'natoms': self.size_n_adj.value,
                                                     'structure': struct,
                                                     'rounding': rounding
                                                      }
        self.makeinfo()

    def clearatoms(self):
        self.atoms = None
        self.pybut.python = None

    def get_atomic_volume(self):
        s = self.list_of_structures[self.structure.get_active()]
        a = self.lattice_const_a.value
        c = self.lattice_const_c.value
        if s == 'fcc':
            return a**3 / 4
        elif s == 'bcc':
            return a**3 / 2
        elif s == 'sc':
            return a**3
        elif s == 'hcp':
            return np.sqrt(3.0)/2 * a * a * c / 2
        elif s == 'graphite':
            return np.sqrt(3.0)/2 * a * a * c / 4
        else:
            raise RuntimeError("Unknown structure: "+s)

    def makeinfo(self):
        """Fill in information field about the atoms.

        Also turns the Wulff construction buttons [Larger] and
        [Smaller] on and off.
        """
        if self.atoms is None:
            self.natoms_label.set_label("-")
            self.dia1_label.set_label("-")
            self.smaller_button.set_sensitive(False)
            self.larger_button.set_sensitive(False)
        else:
            self.natoms_label.set_label(str(len(self.atoms)))
            at_vol = self.get_atomic_volume()
            dia = 2 * (3 * len(self.atoms) * at_vol / (4 * np.pi))**(1.0/3.0)
            self.dia1_label.set_label(_(u"%.1f Å") % (dia,))
            self.smaller_button.set_sensitive(True)
            self.larger_button.set_sensitive(True)
            
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
            
