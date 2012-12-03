"Base class for simulation windows"

import gtk
from gettext import gettext as _
from ase.gui.widgets import oops, pack, help
from ase import Atoms
from ase.constraints import FixAtoms

class Simulation(gtk.Window):
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.gui = gui

    def packtext(self, vbox, text, label=None):
        "Pack an text frame into the window."
        pack(vbox, gtk.Label(""))
        txtframe = gtk.Frame(label)
        txtlbl = gtk.Label(text)
        txtframe.add(txtlbl)
        txtlbl.show()
        pack(vbox, txtframe)
        pack(vbox, gtk.Label(""))

    def packimageselection(self, outerbox, txt1=_(" (rerun simulation)"),
                           txt2=_(" (continue simulation)")):
        "Make the frame for selecting starting config if more than one."
        self.startframe = gtk.Frame(_("Select starting configuration:"))
        pack(outerbox, [self.startframe])
        vbox = gtk.VBox()
        self.startframe.add(vbox)
        vbox.show()
        self.numconfig_format = _("There are currently %i "
                                  "configurations loaded.")
        self.numconfig_label = gtk.Label("")
        pack(vbox, [self.numconfig_label])
        lbl = gtk.Label(_("Choose which one to use as the "
                          "initial configuration"))
        pack(vbox, [lbl])
        self.start_radio_first = gtk.RadioButton(
            None, _("The first configuration %s.") % txt1)
        pack(vbox, [self.start_radio_first])
        self.start_radio_nth = gtk.RadioButton(self.start_radio_first,
                                               _("Configuration number "))
        self.start_nth_adj = gtk.Adjustment(0, 0, 1, 1)
        self.start_nth_spin = gtk.SpinButton(self.start_nth_adj, 0, 0)
        self.start_nth_spin.set_sensitive(False)
        pack(vbox, [self.start_radio_nth, self.start_nth_spin])
        self.start_radio_last = gtk.RadioButton(self.start_radio_first,
            _("The last configuration %s.") % txt2)
        self.start_radio_last.set_active(True)
        pack(vbox, self.start_radio_last)
        self.start_radio_nth.connect("toggled", self.start_radio_nth_toggled)
        self.setupimageselection()
        
    def start_radio_nth_toggled(self, widget):
        self.start_nth_spin.set_sensitive(self.start_radio_nth.get_active())

    def setupimageselection(self):
        "Decide if the start image selection frame should be shown."
        n = self.gui.images.nimages
        if n <= 1:
            self.startframe.hide()
        else:
            self.startframe.show()
            if self.start_nth_adj.value >= n:
                self.start_nth_adj.value = n-1
            self.start_nth_adj.upper = n-1
            self.numconfig_label.set_text(self.numconfig_format % (n,))

    def getimagenumber(self):
        "Get the image number selected in the start image frame."
        nmax = self.gui.images.nimages
        if nmax <= 1:
            return 0
        elif self.start_radio_first.get_active():
            return 0
        elif self.start_radio_nth.get_active():
            return self.start_nth_adj.value
        else:
            assert self.start_radio_last.get_active()
            return nmax-1

    def makebutbox(self, vbox, helptext=None):
        self.buttons = gtk.HButtonBox()
        runbut = gtk.Button(_("Run"))
        runbut.connect('clicked', self.run)
        closebut = gtk.Button(stock=gtk.STOCK_CLOSE)
        closebut.connect('clicked', lambda x: self.destroy())
        for w in (runbut, closebut):
            self.buttons.pack_start(w, 0, 0)
            w.show()
        if helptext:
            helpbut = [help(helptext)]
        else:
            helpbut = []
        pack(vbox, helpbut + [self.buttons], end=True, bottom=True)

    def setup_atoms(self):
        self.atoms = self.get_atoms()
        if self.atoms is None:
            return False
        try:
            self.calculator = self.gui.simulation['calc']
        except KeyError:
            oops(_("No calculator: Use Calculate/Set Calculator on the menu."))
            return False
        self.atoms.set_calculator(self.calculator())
        return True
    
    def get_atoms(self):
        "Make an atoms object from the active image"
        images = self.gui.images
        if images.natoms < 1:
            oops(_("No atoms present"))
            return None
        n = self.getimagenumber()
        natoms = len(images.P[n]) / images.repeat.prod()
        constraint = None
        if not images.dynamic.all():
            constraint = FixAtoms(mask=1-images.dynamic)
        return Atoms(positions=images.P[n,:natoms],
                     symbols=images.Z[:natoms],
                     cell=images.A[n],
                     magmoms=images.M[n,:natoms],
                     tags=images.T[n,:natoms],
                     pbc=images.pbc,
                     constraint=constraint)

    def begin(self, **kwargs):
        if self.gui.simulation.has_key('progress'):
            self.gui.simulation['progress'].begin(**kwargs)

    def end(self):
        if self.gui.simulation.has_key('progress'):
            self.gui.simulation['progress'].end()

    def prepare_store_atoms(self):
        "Informs the gui that the next configuration should be the first."
        self.gui.prepare_new_atoms()
        self.count_steps = 0
        
    def store_atoms(self):
        "Observes the minimization and stores the atoms in the gui."
        self.gui.append_atoms(self.atoms)
        self.count_steps += 1

