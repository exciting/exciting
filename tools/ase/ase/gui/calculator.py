# encoding: utf-8
"calculator.py - module for choosing a calculator."

import gtk
from gettext import gettext as _
import os
import numpy as np
from copy import copy
from ase.gui.setupwindow import SetupWindow
from ase.gui.progress import DefaultProgressIndicator, GpawProgressIndicator
from ase.gui.widgets import pack, oops, cancel_apply_ok
from ase import Atoms
from ase.data import chemical_symbols
import ase

# Asap and GPAW may be imported if selected.

introtext = _("""\
To make most calculations on the atoms, a Calculator object must first
be associated with it.  ASE supports a number of calculators, supporting
different elements, and implementing different physical models for the
interatomic interactions.\
""")

# Informational text about the calculators
lj_info_txt = _("""\
The Lennard-Jones pair potential is one of the simplest
possible models for interatomic interactions, mostly
suitable for noble gasses and model systems.

Interactions are described by an interaction length and an
interaction strength.\
""")

emt_info_txt = _("""\
The EMT potential is a many-body potential, giving a
good description of the late transition metals crystalling
in the FCC crystal structure.  The elements described by the
main set of EMT parameters are Al, Ni, Cu, Pd, Ag, Pt, and
Au, the Al potential is however not suitable for materials
science application, as the stacking fault energy is wrong.

A number of parameter sets are provided.

<b>Default parameters:</b>

The default EMT parameters, as published in K. W. Jacobsen,
P. Stoltze and J. K. Nørskov, <i>Surf. Sci.</i> <b>366</b>, 394 (1996).

<b>Alternative Cu, Ag and Au:</b>

An alternative set of parameters for Cu, Ag and Au,
reoptimized to experimental data including the stacking
fault energies by Torben Rasmussen (partly unpublished).

<b>Ruthenium:</b>

Parameters for Ruthenium, as published in J. Gavnholt and
J. Schiøtz, <i>Phys. Rev. B</i> <b>77</b>, 035404 (2008).

<b>Metallic glasses:</b>

Parameters for MgCu and CuZr metallic glasses. MgCu
parameters are in N. P. Bailey, J. Schiøtz and
K. W. Jacobsen, <i>Phys. Rev. B</i> <b>69</b>, 144205 (2004).
CuZr in A. Paduraru, A. Kenoufi, N. P. Bailey and
J. Schiøtz, <i>Adv. Eng. Mater.</i> <b>9</b>, 505 (2007).
""")

aseemt_info_txt = _("""\
The EMT potential is a many-body potential, giving a
good description of the late transition metals crystalling
in the FCC crystal structure.  The elements described by the
main set of EMT parameters are Al, Ni, Cu, Pd, Ag, Pt, and
Au.  In addition, this implementation allows for the use of
H, N, O and C adatoms, although the description of these is
most likely not very good.

<b>This is the ASE implementation of EMT.</b> For large
simulations the ASAP implementation is more suitable; this
implementation is mainly to make EMT available when ASAP is
not installed.
""")

brenner_info_txt = _("""\
The Brenner potential is a reactive bond-order potential for
carbon and hydrocarbons.  As a bond-order potential, it takes
into account that carbon orbitals can hybridize in different
ways, and that carbon can form single, double and triple
bonds.  That the potential is reactive means that it can
handle gradual changes in the bond order as chemical bonds
are formed or broken.

The Brenner potential is implemented in Asap, based on a
C implentation published at http://www.rahul.net/pcm/brenner/ .

The potential is documented here:
  Donald W Brenner, Olga A Shenderova, Judith A Harrison,
  Steven J Stuart, Boris Ni and Susan B Sinnott:
  "A second-generation reactive empirical bond order (REBO)
  potential energy expression for hydrocarbons",
  J. Phys.: Condens. Matter 14 (2002) 783-802.
  doi: 10.1088/0953-8984/14/4/312
""")


gpaw_info_txt = _("""\
GPAW implements Density Functional Theory using a
<b>G</b>rid-based real-space representation of the wave
functions, and the <b>P</b>rojector <b>A</b>ugmented <b>W</b>ave
method for handling the core regions.  
""")

aims_info_txt = _("""\
FHI-aims is an external package implementing density 
functional theory and quantum chemical methods using 
all-electron methods and a numeric local orbital basis set. 
For full details, see http://www.fhi-berlin.mpg.de/aims/ 
or Comp. Phys. Comm. v180 2175 (2009). The ASE 
documentation contains information on the keywords and 
functionalities available within this interface. 
""")

aims_pbc_warning_text = _("""\
WARNING:
Your system seems to have more than zero but less than 
three periodic dimensions. Please check that this is 
really what you want to compute. Assuming full 
3D periodicity for this calculator.""")

vasp_info_txt = _("""\
VASP is an external package implementing density 
functional functional theory using pseudopotentials 
or the projector-augmented wave method together 
with a plane wave basis set. For full details, see
http://cms.mpi.univie.ac.at/vasp/vasp/
""")

emt_parameters = (
    (_("Default (Al, Ni, Cu, Pd, Ag, Pt, Au)"), None),
    (_("Alternative Cu, Ag and Au"), "EMTRasmussenParameters"),
    (_("Ruthenium"), "EMThcpParameters"),
    (_("CuMg and CuZr metallic glass"), "EMTMetalGlassParameters")
    )

class SetCalculator(SetupWindow):
    "Window for selecting a calculator."

    # List the names of the radio button attributes
    radios = ("none", "lj", "emt", "aseemt", "brenner", "gpaw", "aims", "vasp")
    # List the names of the parameter dictionaries
    paramdicts = ("lj_parameters","gpaw_parameters","aims_parameters",)
    # The name used to store parameters on the gui object
    classname = "SetCalculator"
    
    def __init__(self, gui):
        SetupWindow.__init__(self)
        self.set_title(_("Select calculator"))
        vbox = gtk.VBox()
        
        # Intoductory text
        self.packtext(vbox, introtext)
        
        pack(vbox, [gtk.Label(_("Calculator:"))])

        # No calculator (the default)
        self.none_radio = gtk.RadioButton(None, _("None"))
        pack(vbox, [self.none_radio])

        # Lennard-Jones
        self.lj_radio = gtk.RadioButton(self.none_radio,
                                        _("Lennard-Jones (ASAP)"))
        self.lj_setup = gtk.Button(_("Setup"))
        self.lj_info = InfoButton(lj_info_txt)
        self.lj_setup.connect("clicked", self.lj_setup_window)
        self.pack_line(vbox, self.lj_radio, self.lj_setup, self.lj_info)

        # EMT
        self.emt_radio = gtk.RadioButton(
            self.none_radio, _("EMT - Effective Medium Theory (ASAP)"))
        self.emt_setup = gtk.combo_box_new_text()
        self.emt_param_info = {}
        for p in emt_parameters:
            self.emt_setup.append_text(p[0])
            self.emt_param_info[p[0]] = p[1]
        self.emt_setup.set_active(0)
        self.emt_info = InfoButton(emt_info_txt)
        self.pack_line(vbox, self.emt_radio, self.emt_setup, self.emt_info)

        # EMT (ASE implementation)
        self.aseemt_radio = gtk.RadioButton(
            self.none_radio, _("EMT - Effective Medium Theory (ASE)"))
        self.aseemt_info = InfoButton(aseemt_info_txt)
        self.pack_line(vbox, self.aseemt_radio, None, self.aseemt_info)

        # Brenner potential
        self.brenner_radio = gtk.RadioButton(
            self.none_radio, _("Brenner Potential (ASAP)"))
        self.brenner_info = InfoButton(brenner_info_txt)
        self.pack_line(vbox, self.brenner_radio, None, self.brenner_info)
        
        # GPAW
        self.gpaw_radio = gtk.RadioButton(self.none_radio,
                                          _("Density Functional Theory (GPAW)")
                                          )
        self.gpaw_setup = gtk.Button(_("Setup"))
        self.gpaw_info = InfoButton(gpaw_info_txt)
        self.gpaw_setup.connect("clicked", self.gpaw_setup_window)
        self.pack_line(vbox, self.gpaw_radio, self.gpaw_setup, self.gpaw_info)
        
        # FHI-aims
        self.aims_radio = gtk.RadioButton(self.none_radio, 
                                          _("Density Functional Theory "
                                            "(FHI-aims)"))
        self.aims_setup = gtk.Button(_("Setup"))
        self.aims_info = InfoButton(aims_info_txt)
        self.aims_setup.connect("clicked", self.aims_setup_window)
        self.pack_line(vbox, self.aims_radio, self.aims_setup, self.aims_info)
        
        # VASP
        self.vasp_radio = gtk.RadioButton(self.none_radio, 
                                          _("Density Functional Theory "
                                            "(VASP)"))
        self.vasp_setup = gtk.Button(_("Setup"))
        self.vasp_info = InfoButton(vasp_info_txt)
        self.vasp_setup.connect("clicked", self.vasp_setup_window)
        self.pack_line(vbox, self.vasp_radio, self.vasp_setup, self.vasp_info)

        # Buttons etc.
        pack(vbox, gtk.Label(""))
        buts = cancel_apply_ok(cancel=lambda widget: self.destroy(),
                               apply=self.apply,
                               ok=self.ok)
        pack(vbox, [buts], end=True, bottom=True)
        self.check = gtk.CheckButton(_("Check that the calculator is "
                                       "reasonable."))
        self.check.set_active(True)
        fr = gtk.Frame()
        fr.add(self.check)
        fr.show_all()
        pack(vbox, [fr], end=True, bottom=True)
        
        # Finalize setup
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui
        self.load_state()
        
    def pack_line(self, box, radio, setup, info):
        hbox = gtk.HBox()
        hbox.pack_start(radio, 0, 0)
        hbox.pack_start(gtk.Label("  "), 0, 0)
        hbox.pack_end(info, 0, 0)
        if setup is not None:
            radio.connect("toggled", self.radio_toggled, setup)
            setup.set_sensitive(False)
            hbox.pack_end(setup, 0, 0)
        hbox.show_all()
        box.pack_start(hbox, 0, 0)

    def radio_toggled(self, radio, button):
        button.set_sensitive(radio.get_active())

    def lj_setup_window(self, widget):
        if not self.get_atoms():
            return
        lj_param = getattr(self, "lj_parameters", None)
        LJ_Window(self, lj_param, "lj_parameters")
        # When control is retuned, self.lj_parameters has been set.
        
    def gpaw_setup_window(self, widget):
        if not self.get_atoms():
            return
        gpaw_param = getattr(self, "gpaw_parameters", None)
        GPAW_Window(self, gpaw_param, "gpaw_parameters")
        # When control is retuned, self.gpaw_parameters has been set.
        
    def aims_setup_window(self, widget):
        if not self.get_atoms():
            return
        aims_param = getattr(self, "aims_parameters", None)
        AIMS_Window(self, aims_param, "aims_parameters")        
        # When control is retuned, self.aims_parameters has been set.

    def vasp_setup_window(self, widget):
        if not self.get_atoms():
            return
        vasp_param = getattr(self, "vasp_parameters", None)
        VASP_Window(self, vasp_param, "vasp_parameters")
        # When control is retuned, self.vasp_parameters has been set.

    def get_atoms(self):
        "Make an atoms object from the active frame"
        images = self.gui.images
        frame = self.gui.frame
        if images.natoms < 1:
            oops(_("No atoms present"))
            return False
        self.atoms = Atoms(positions=images.P[frame],
                           symbols=images.Z,
                           cell=images.A[frame],
                           pbc=images.pbc,
                           magmoms=images.M[frame])
        if not images.dynamic.all(): 
            from ase.constraints import FixAtoms
            self.atoms.set_constraint(FixAtoms(mask=1-images.dynamic))
        return True

    def apply(self, *widget):
        if self.do_apply():
            self.save_state()
            return True
        else:
            return False
        
    def do_apply(self):
        nochk = not self.check.get_active()
        self.gui.simulation["progress"] = DefaultProgressIndicator()
        if self.none_radio.get_active():
            self.gui.simulation['calc'] = None
            return True
        elif self.lj_radio.get_active():
            if nochk or self.lj_check():
                self.choose_lj()
                return True
        elif self.emt_radio.get_active():
            if nochk or self.emt_check():
                self.choose_emt()
                return True
        elif self.aseemt_radio.get_active():
            if nochk or self.aseemt_check():
                self.choose_aseemt()
                return True
        elif self.brenner_radio.get_active():
            if nochk or self.brenner_check():
                self.choose_brenner()
                return True
        elif self.gpaw_radio.get_active():
            if nochk or self.gpaw_check():
                self.choose_gpaw()
                return True
        elif self.aims_radio.get_active():
            if nochk or self.aims_check():
                self.choose_aims()
                return True
        elif self.vasp_radio.get_active():
            if nochk or self.vasp_check():
                self.choose_vasp()
                return True  
        return False

    def ok(self, *widget):
        if self.apply():
            self.destroy()

    def save_state(self):
        state = {}
        for r in self.radios:
            radiobutton = getattr(self, r+"_radio")
            if radiobutton.get_active():
                state["radio"] = r
        state["emtsetup"] = self.emt_setup.get_active()
        state["check"] = self.check.get_active()
        for p in self.paramdicts:
            if hasattr(self, p):
                state[p] = getattr(self, p)
        self.gui.module_state[self.classname] = state

    def load_state(self):
        try:
            state = self.gui.module_state[self.classname]
        except KeyError:
            return
        r = state["radio"]
        radiobutton = getattr(self, r + "_radio")
        radiobutton.set_active(True)
        self.emt_setup.set_active(state["emtsetup"])
        self.check.set_active(state["check"])
        for p in self.paramdicts:
            if state.has_key(p):
                setattr(self, p, state[p])
            
    def lj_check(self):
        try:
            import asap3
        except ImportError:
            oops(_("ASAP is not installed. (Failed to import asap3)"))
            return False
        if not hasattr(self, "lj_parameters"):
            oops(_("You must set up the Lennard-Jones parameters"))
            return False
        try:
            self.atoms.set_calculator(asap3.LennardJones(**self.lj_parameters))
        except (asap3.AsapError, TypeError, ValueError), e:
            oops(_("Could not create useful Lennard-Jones calculator."),
                 str(e))
            return False
        return True

    def choose_lj(self):
        # Define a function on the fly!
        import asap3
        def lj_factory(p=self.lj_parameters, lj=asap3.LennardJones):
            return lj(**p)
        self.gui.simulation["calc"] = lj_factory

    def emt_get(self):
        import asap3
        provider_name = self.emt_setup.get_active_text()
        provider =  self.emt_param_info[provider_name]
        if provider is not None:
            provider = getattr(asap3, provider)
        return (asap3.EMT, provider, asap3)
                                      
    def emt_check(self):
        if not self.get_atoms():
            return False
        try:
            emt, provider, asap3 = self.emt_get()
        except ImportError:
            oops(_("ASAP is not installed. (Failed to import asap3)"))
            return False
        try:
            if provider is not None:
                self.atoms.set_calculator(emt(provider()))
            else:
                self.atoms.set_calculator(emt())
        except (asap3.AsapError, TypeError, ValueError), e:
            oops(_("Could not attach EMT calculator to the atoms."),
                 str(e))
            return False
        return True

    def choose_emt(self):
        emt, provider, asap3 = self.emt_get()
        if provider is None:
            emt_factory = emt
        else:
            def emt_factory(emt=emt, prov=provider):
                return emt(prov())
        self.gui.simulation["calc"] = emt_factory

    def aseemt_check(self):
        return self.element_check("ASE EMT", ['H', 'Al', 'Cu', 'Ag', 'Au',
                                              'Ni', 'Pd', 'Pt', 'C', 'N', 'O'])

    def brenner_check(self):
        try:
            import asap3
        except ImportError:
            oops(_("ASAP is not installed. (Failed to import asap3)"))
            return False
        return self.element_check("Brenner potential", ['H', 'C', 'Si'])

    def choose_brenner(self):
        import asap3
        self.gui.simulation["calc"] = asap3.BrennerPotential

    def choose_aseemt(self):
        import ase.calculators.emt 
        self.gui.simulation["calc"] = ase.calculators.emt.EMT
        # In case Asap has been imported
        ase.calculators.emt.EMT.disabled = False

    def gpaw_check(self):
        try:
            import gpaw
        except ImportError:
            oops(_("GPAW is not installed. (Failed to import gpaw)"))
            return False
        if not hasattr(self, "gpaw_parameters"):
            oops(_("You must set up the GPAW parameters"))
            return False
        return True

    def choose_gpaw(self):
        # This reuses the same GPAW object.
        try:
            import gpaw
        except ImportError:
            oops(_("GPAW is not installed. (Failed to import gpaw)"))
            return False
        p = self.gpaw_parameters
        use = ["xc", "kpts", "mode"]
        if p["use_h"]:
            use.append("h")
        else:
            use.append("gpts")
        if p["mode"] == "lcao":
            use.append("basis")
        gpaw_param = {}
        for s in use:
            gpaw_param[s] = p[s]
        if p["use mixer"]:
            mx = getattr(gpaw, p["mixer"])
            mx_args = {}
            mx_arg_n = ["beta", "nmaxold", "weight"]
            if p["mixer"] == "MixerDiff":
                mx_arg_n.extend(["beta_m", "nmaxold_m", "weight_m"])
            for s in mx_arg_n:
                mx_args[s] = p[s]
            gpaw_param["mixer"] = mx(**mx_args)
        progress = GpawProgressIndicator()
        self.gui.simulation["progress"] = progress
        gpaw_param["txt"] = progress.get_gpaw_stream()
        gpaw_calc = gpaw.GPAW(**gpaw_param)
        def gpaw_factory(calc = gpaw_calc):
            return calc
        self.gui.simulation["calc"] = gpaw_factory
                
    def aims_check(self):
        if not hasattr(self, "aims_parameters"):
            oops(_("You must set up the FHI-aims parameters"))
            return False
        return True

    def choose_aims(self):
        param = self.aims_parameters
        from ase.calculators.aims import Aims
        calc_aims = Aims(**param)
        def aims_factory(calc = calc_aims):
            return calc
        self.gui.simulation["calc"] = aims_factory

    def vasp_check(self):
        if not hasattr(self, "vasp_parameters"):
            oops(_("You must set up the VASP parameters"))
            return False
        return True

    def choose_vasp(self):
        param = self.vasp_parameters
        from ase.calculators.vasp import Vasp
        calc_vasp = Vasp(**param)
        def vasp_factory(calc = calc_vasp):
            return calc
        self.gui.simulation["calc"] = vasp_factory

    def element_check(self, name, elements):
        "Check that all atoms are allowed"
        elements = [ase.data.atomic_numbers[s] for s in elements]
        elements_dict = {}
        for e in elements:
            elements_dict[e] = True
        if not self.get_atoms():
            return False
        try:
            for e in self.atoms.get_atomic_numbers():
                elements_dict[e]
        except KeyError:
            oops(_("Element %(sym)s not allowed by the '%(name)s' calculator")
                 % dict(sym=ase.data.chemical_symbols[e], name=name))
            return False
        return True
    
class InfoButton(gtk.Button):
    def __init__(self, txt):
        gtk.Button.__init__(self, _("Info"))
        self.txt = txt
        self.connect('clicked', self.run)

    def run(self, widget):
        dialog = gtk.MessageDialog(flags=gtk.DIALOG_MODAL,
                                   type=gtk.MESSAGE_INFO,
                                   buttons=gtk.BUTTONS_CLOSE)
        dialog.set_markup(self.txt)
        dialog.connect('response', lambda x, y: dialog.destroy())
        dialog.show()


class LJ_Window(gtk.Window):
    def __init__(self, owner, param, attrname):
        gtk.Window.__init__(self)
        self.set_title(_("Lennard-Jones parameters"))
        self.owner = owner
        self.attrname = attrname
        atoms = owner.atoms
        atnos = atoms.get_atomic_numbers()
        found = {}
        for z in atnos:
            found[z] = True
        self.present = found.keys()
        self.present.sort()  # Sorted list of atomic numbers
        nelem = len(self.present)
        vbox = gtk.VBox()
        label = gtk.Label(_("Specify the Lennard-Jones parameters here"))
        pack(vbox, [label])
        pack(vbox, gtk.Label(""))
        pack(vbox, [gtk.Label(_("Epsilon (eV):"))])
        tbl, self.epsilon_adj = self.makematrix(self.present)
        pack(vbox, [tbl])
        pack(vbox, gtk.Label(""))
        pack(vbox, [gtk.Label(_(u"Sigma (Å):"))])
        tbl, self.sigma_adj = self.makematrix(self.present)
        pack(vbox, [tbl])
        # TRANSLATORS: Shift roughly means adjust (about a potential)
        self.modif = gtk.CheckButton(_("Shift to make smooth at cutoff"))
        self.modif.set_active(True)
        pack(vbox, gtk.Label(""))
        pack(vbox, self.modif)
        pack(vbox, gtk.Label(""))
        butbox = gtk.HButtonBox()
        cancel_but = gtk.Button(stock=gtk.STOCK_CANCEL)
        cancel_but.connect('clicked', lambda widget: self.destroy())
        ok_but = gtk.Button(stock=gtk.STOCK_OK)
        ok_but.connect('clicked', self.ok)
        butbox.pack_start(cancel_but, 0, 0)
        butbox.pack_start(ok_but, 0, 0)
        butbox.show_all()
        pack(vbox, [butbox], end=True, bottom=True)
        vbox.show()
        self.add(vbox)

        # Now, set the parameters
        if param and param['elements'] == self.present:
            self.set_param(self.epsilon_adj, param["epsilon"], nelem)
            self.set_param(self.sigma_adj, param["sigma"], nelem)
            self.modif.set_active(param["modified"])

        self.show()
        self.grab_add()  # Lock all other windows
        
    def makematrix(self, present):
        nelem = len(present)
        adjdict = {}
        tbl = gtk.Table(2+nelem, 2+nelem)
        for i in range(nelem):
            s = chemical_symbols[present[i]]
            tbl.attach(gtk.Label(" " + str(present[i])), 0, 1, i, i+1)
            tbl.attach(gtk.Label("  "+s+" "), 1, 2, i, i+1)
            tbl.attach(gtk.Label(str(present[i])), i+2, i+3, 1+nelem, 2+nelem)
            tbl.attach(gtk.Label(s), i+2, i+3, nelem, 1+nelem)
            for j in range(i+1):
                adj = gtk.Adjustment(1.0, 0.0, 100.0, 0.1)
                spin = gtk.SpinButton(adj, 0.1, 3)
                tbl.attach(spin, 2+j, 3+j, i, i+1)
                adjdict[(i,j)] = adj
        tbl.show_all()
        return tbl, adjdict
    
    def set_param(self, adj, params, n):
        for i in range(n):
            for j in range(n):
                if j <= i:
                    adj[(i,j)].value = params[i,j]

    def get_param(self, adj, params, n):
        for i in range(n):
            for j in range(n):
                if j <= i:
                    params[i,j] = params[j,i] = adj[(i,j)].value


    def destroy(self):
        self.grab_remove()
        gtk.Window.destroy(self)

    def ok(self, *args):
        params = {}
        params["elements"] = copy(self.present)
        n = len(self.present)
        eps = np.zeros((n,n))
        self.get_param(self.epsilon_adj, eps, n)
        sigma = np.zeros((n,n))
        self.get_param(self.sigma_adj, sigma, n)
        params["epsilon"] = eps
        params["sigma"] = sigma
        params["modified"] = self.modif.get_active()
        setattr(self.owner, self.attrname, params)
        self.destroy()


class GPAW_Window(gtk.Window):
    gpaw_xc_list = ['LDA', 'PBE', 'RPBE', 'revPBE']
    gpaw_xc_default = 'PBE'
    def __init__(self, owner, param, attrname):
        gtk.Window.__init__(self)
        self.set_title(_("GPAW parameters"))
        self.owner = owner
        self.attrname = attrname
        atoms = owner.atoms
        self.ucell = atoms.get_cell()
        self.size = tuple([self.ucell[i,i] for i in range(3)])
        self.pbc = atoms.get_pbc()
        self.orthogonal = self.isorthogonal(self.ucell)
        self.natoms = len(atoms)
        
        vbox = gtk.VBox()
        #label = gtk.Label("Specify the GPAW parameters here")
        #pack(vbox, [label])

        # Print some info
        txt = _("%i atoms.\n") % (self.natoms,)
        if self.orthogonal:
            txt += _(u"Orthogonal unit cell: %.2f x %.2f x %.2f Å.") % self.size
        else:
            txt += _("Non-orthogonal unit cell:\n")
            txt += str(self.ucell)
        pack(vbox, [gtk.Label(txt)])
        
        # XC potential
        self.xc = gtk.combo_box_new_text()
        for i, x in enumerate(self.gpaw_xc_list):
            self.xc.append_text(x)
            if x == self.gpaw_xc_default:
                self.xc.set_active(i)
        pack(vbox, [gtk.Label(_("Exchange-correlation functional: ")),
                    self.xc])
        
        # Grid spacing
        self.radio_h = gtk.RadioButton(None, _("Grid spacing"))
        self.h = gtk.Adjustment(0.18, 0.0, 1.0, 0.01)
        self.h_spin = gtk.SpinButton(self.h, 0, 2)
        pack(vbox, [self.radio_h, gtk.Label(" h = "), self.h_spin,
                    gtk.Label(_(u"Å"))])
        self.radio_gpts = gtk.RadioButton(self.radio_h, _("Grid points"))
        self.gpts = []
        self.gpts_spin = []
        for i in range(3):
            g = gtk.Adjustment(4, 4, 1000, 4)
            s = gtk.SpinButton(g, 0, 0)
            self.gpts.append(g)
            self.gpts_spin.append(s)
        self.gpts_hlabel = gtk.Label("")
        self.gpts_hlabel_format = _(u"h<sub>eff</sub> = (%.3f, %.3f, %.3f) Å")
        pack(vbox, [self.radio_gpts, gtk.Label(" gpts = ("), self.gpts_spin[0],
                    gtk.Label(", "), self.gpts_spin[1], gtk.Label(", "),
                    self.gpts_spin[2], gtk.Label(")  "), self.gpts_hlabel])
        self.radio_h.connect("toggled", self.radio_grid_toggled)
        self.radio_gpts.connect("toggled", self.radio_grid_toggled)
        self.radio_grid_toggled(None)
        for g in self.gpts:
            g.connect("value-changed", self.gpts_changed)
        self.h.connect("value-changed", self.h_changed)
        
        # K-points
        self.kpts = []
        self.kpts_spin = []
        for i in range(3):
            if self.pbc[i] and self.orthogonal:
                default = np.ceil(20.0 / self.size[i])
            else:
                default = 1
            g = gtk.Adjustment(default, 1, 100, 1)
            s = gtk.SpinButton(g, 0, 0)
            self.kpts.append(g)
            self.kpts_spin.append(s)
            if not self.pbc[i]:
                s.set_sensitive(False)
            g.connect("value-changed", self.k_changed)
        pack(vbox, [gtk.Label(_("k-points  k = (")), self.kpts_spin[0],
                    gtk.Label(", "), self.kpts_spin[1], gtk.Label(", "),
                    self.kpts_spin[2], gtk.Label(")")])
        self.kpts_label = gtk.Label("")
        self.kpts_label_format = _(u"k-points x size:  (%.1f, %.1f, %.1f) Å")
        pack(vbox, [self.kpts_label])
        self.k_changed()
        
        # Spin polarized
        self.spinpol = gtk.CheckButton(_("Spin polarized"))
        pack(vbox, [self.spinpol])
        pack(vbox, gtk.Label(""))

        # Mode and basis functions
        self.mode = gtk.combo_box_new_text()
        self.mode.append_text(_("FD - Finite Difference (grid) mode"))
        self.mode.append_text(_("LCAO - Linear Combination of Atomic "
                                "Orbitals"))
        self.mode.set_active(0)
        pack(vbox, [gtk.Label(_("Mode: ")), self.mode])
        self.basis = gtk.combo_box_new_text()
        self.basis.append_text(_("sz - Single Zeta"))
        self.basis.append_text(_("szp - Single Zeta polarized"))
        self.basis.append_text(_("dzp - Double Zeta polarized"))
        self.basis.set_active(2) # dzp
        pack(vbox, [gtk.Label(_("Basis functions: ")), self.basis])
        pack(vbox, gtk.Label(""))
        self.mode.connect("changed", self.mode_changed)
        self.mode_changed()
        
        # Mixer
        self.use_mixer = gtk.CheckButton(_("Non-standard mixer parameters"))
        pack(vbox, [self.use_mixer])
        self.radio_mixer = gtk.RadioButton(None, "Mixer   ")
        self.radio_mixersum = gtk.RadioButton(self.radio_mixer, "MixerSum   ")
        self.radio_mixerdiff = gtk.RadioButton(self.radio_mixer, "MixerDiff")
        pack(vbox, [self.radio_mixer, self.radio_mixersum,
                    self.radio_mixerdiff])
        self.beta_adj = gtk.Adjustment(0.25, 0.0, 1.0, 0.05)
        self.beta_spin = gtk.SpinButton(self.beta_adj, 0, 2)
        self.nmaxold_adj = gtk.Adjustment(3, 1, 10, 1)
        self.nmaxold_spin = gtk.SpinButton(self.nmaxold_adj, 0, 0)
        self.weight_adj = gtk.Adjustment(50, 1, 500, 1)
        self.weight_spin = gtk.SpinButton(self.weight_adj, 0, 0)
        pack(vbox, [gtk.Label("beta = "), self.beta_spin,
                    gtk.Label("  nmaxold = "), self.nmaxold_spin,
                    gtk.Label("  weight = "), self.weight_spin])
        self.beta_m_adj = gtk.Adjustment(0.70, 0.0, 1.0, 0.05)
        self.beta_m_spin = gtk.SpinButton(self.beta_m_adj, 0, 2)
        self.nmaxold_m_adj = gtk.Adjustment(2, 1, 10, 1)
        self.nmaxold_m_spin = gtk.SpinButton(self.nmaxold_m_adj, 0, 0)
        self.weight_m_adj = gtk.Adjustment(10, 1, 500, 1)
        self.weight_m_spin = gtk.SpinButton(self.weight_m_adj, 0, 0)
        pack(vbox, [gtk.Label("beta_m = "), self.beta_m_spin,
                    gtk.Label("  nmaxold_m = "), self.nmaxold_m_spin,
                    gtk.Label("  weight_m = "), self.weight_m_spin])
        for but in (self.spinpol, self.use_mixer, self.radio_mixer,
                    self.radio_mixersum, self.radio_mixerdiff):
            but.connect("clicked", self.mixer_changed)
        self.mixer_changed()
        
        # Eigensolver
        # Poisson-solver
        
        vbox.show()
        self.add(vbox)

        # Buttons at the bottom
        pack(vbox, gtk.Label(""))
        butbox = gtk.HButtonBox()
        cancel_but = gtk.Button(stock=gtk.STOCK_CANCEL)
        cancel_but.connect('clicked', lambda widget: self.destroy())
        ok_but = gtk.Button(stock=gtk.STOCK_OK)
        ok_but.connect('clicked', self.ok)
        butbox.pack_start(cancel_but, 0, 0)
        butbox.pack_start(ok_but, 0, 0)
        butbox.show_all()
        pack(vbox, [butbox], end=True, bottom=True)

        # Set stored parameters
        if param:
            self.xc.set_active(param["xc#"])
            if param["use_h"]:
                self.radio_h.set_active(True)
            else:
                self.radio_gpts.set_active(True)
            for i in range(3):
                self.gpts[i].value = param["gpts"][i]
                self.kpts[i].value = param["kpts"][i]
            self.spinpol.set_active(param["spinpol"])
            self.mode.set_active(param["mode#"])
            self.basis.set_active(param["basis#"])
            self.use_mixer.set_active(param["use mixer"])
            getattr(self, "radio_" + param["mixer"].lower()).set_active(True)
            for t in ("beta", "nmaxold", "weight", "beta_m", "nmaxold_m",
                      "weight_m"):                    
                getattr(self, t+"_adj").value = param[t]

        self.show()
        self.grab_add()  # Lock all other windows

    def radio_grid_toggled(self, widget):
        hmode = self.radio_h.get_active()
        self.h_spin.set_sensitive(hmode)
        for s in self.gpts_spin:
            s.set_sensitive(not hmode)
        self.gpts_changed()

    def gpts_changed(self, *args):
        if self.radio_gpts.get_active():
            g = np.array([int(g.value) for g in self.gpts])
            size = np.array([self.ucell[i,i] for i in range(3)])
            txt = self.gpts_hlabel_format % tuple(size / g)
            self.gpts_hlabel.set_markup(txt)
        else:
            self.gpts_hlabel.set_markup("")

    def h_changed(self, *args):
        h = self.h.value
        for i in range(3):
            g = 4 * round(self.ucell[i,i] / (4*h))
            self.gpts[i].value = g

    def k_changed(self, *args):
        size = [self.kpts[i].value * np.sqrt(np.vdot(self.ucell[i],
                                                     self.ucell[i]))
                for i in range(3)]
        self.kpts_label.set_text(self.kpts_label_format % tuple(size))

    def mode_changed(self, *args):
        self.basis.set_sensitive(self.mode.get_active() == 1)

    def mixer_changed(self, *args):
        radios = (self.radio_mixer, self.radio_mixersum, self.radio_mixerdiff)
        spin1 = (self.beta_spin, self.nmaxold_spin, self.weight_spin)
        spin2 = (self.beta_m_spin, self.nmaxold_m_spin, self.weight_m_spin)
        if self.use_mixer.get_active():
            # Mixer parameters can be specified.
            if self.spinpol.get_active():
                self.radio_mixer.set_sensitive(False)
                self.radio_mixersum.set_sensitive(True)
                self.radio_mixerdiff.set_sensitive(True)
                if self.radio_mixer.get_active():
                    self.radio_mixersum.set_active(True)
            else:
                self.radio_mixer.set_sensitive(True)
                self.radio_mixersum.set_sensitive(False)
                self.radio_mixerdiff.set_sensitive(False)
                self.radio_mixer.set_active(True)
            if self.radio_mixerdiff.get_active():
                active = spin1 + spin2
                passive = ()
            else:
                active = spin1
                passive = spin2
            for widget in active:
                widget.set_sensitive(True)
            for widget in passive:
                widget.set_sensitive(False)
        else:
            # No mixer parameters
            for widget in radios + spin1 + spin2:
                widget.set_sensitive(False)
                
    def isorthogonal(self, matrix):
        ortho = True
        for i in range(3):
            for j in range(3):
                if i != j and matrix[i][j] != 0.0:
                    ortho = False
        return ortho

    def ok(self, *args):
        param = {}
        param["xc"] = self.xc.get_active_text()
        param["xc#"] = self.xc.get_active()
        param["use_h"] = self.radio_h.get_active()
        param["h"] = self.h.value
        param["gpts"] = [int(g.value) for g in self.gpts]
        param["kpts"] = [int(k.value) for k in self.kpts]
        param["spinpol"] = self.spinpol.get_active()
        param["mode"] = self.mode.get_active_text().split()[0].lower()
        param["mode#"] = self.mode.get_active()
        param["basis"] = self.basis.get_active_text().split()[0].lower()
        param["basis#"] = self.basis.get_active()
        param["use mixer"] = self.use_mixer.get_active()
        if self.radio_mixer.get_active():
            m = "Mixer"
        elif self.radio_mixersum.get_active():
            m = "MixerSum"
        else:
            assert self.radio_mixerdiff.get_active()
            m = "MixerDiff"
        param["mixer"] = m
        for t in ("beta", "nmaxold", "weight", "beta_m", "nmaxold_m",
                  "weight_m"):
            param[t] = getattr(self, t+"_adj").value
        setattr(self.owner, self.attrname, param)
        self.destroy()

class AIMS_Window(gtk.Window):
    aims_xc_cluster = ['pw-lda','pz-lda','pbe','pbesol','rpbe','revpbe',
                    'blyp','am05','b3lyp','hse03','hse06','pbe0','pbesol0',
                    'hf','mp2']
    aims_xc_periodic = ['pw-lda','pz-lda','pbe','pbesol','rpbe','revpbe',
                        'blyp','am05']
    aims_xc_default = 'pbe'
    aims_relativity_list = ['none','atomic_zora','zora']
    aims_keyword_gui_list = ['xc','vdw_correction_hirshfeld','k_grid','spin','charge','relativistic',
                             'sc_accuracy_etot','sc_accuracy_eev','sc_accuracy_rho','sc_accuracy_forces',
                             'compute_forces','run_command','species_dir','default_initial_moment']
    def __init__(self, owner, param, attrname):
        self.owner = owner
        self.attrname = attrname
        atoms = owner.atoms
        self.periodic = atoms.get_pbc().all()
        if not self.periodic and atoms.get_pbc().any():
            aims_periodic_warning = True
            self.periodic = True 
        else:
            aims_periodic_warning = False
        from ase.calculators.aims import float_keys,exp_keys,string_keys,int_keys,bool_keys,list_keys,input_keys
        self.aims_keyword_list =float_keys+exp_keys+string_keys+int_keys+bool_keys+list_keys+input_keys
        self.expert_keywords = []

        natoms = len(atoms)
        gtk.Window.__init__(self)
        self.set_title(_("FHI-aims parameters"))
        vbox = gtk.VBox()
        vbox.set_border_width(5)
        # Print some info
        txt = _("%i atoms.\n") % (natoms)
        if self.periodic:
            self.ucell = atoms.get_cell()
            txt += _("Periodic geometry, unit cell is:\n")
            for i in range(3):
                txt += "(%8.3f %8.3f %8.3f)\n" % (self.ucell[i][0], self.ucell[i][1], self.ucell[i][2])
            self.xc_list = self.aims_xc_periodic
        else:
            txt += _("Non-periodic geometry.\n")
            self.xc_list = self.aims_xc_cluster
        pack(vbox, [gtk.Label(txt)])

        # XC functional & dispersion correction
        self.xc = gtk.combo_box_new_text()
        self.xc_setup = False
        self.TS = gtk.CheckButton(_("Hirshfeld-based dispersion correction"))
        pack(vbox, [gtk.Label(_("Exchange-correlation functional: ")),self.xc])
        pack(vbox, [self.TS])
        pack(vbox, [gtk.Label("")])
        
        # k-grid?
        if self.periodic:
            self.kpts = []
            self.kpts_spin = []
            for i in range(3):
                default = np.ceil(20.0 / np.sqrt(np.vdot(self.ucell[i],self.ucell[i])))
                g = gtk.Adjustment(default, 1, 100, 1)
                s = gtk.SpinButton(g, 0, 0)
                self.kpts.append(g)
                self.kpts_spin.append(s)
                g.connect("value-changed", self.k_changed)
            pack(vbox, [gtk.Label(_("k-points  k = (")), self.kpts_spin[0],
                        gtk.Label(", "), self.kpts_spin[1], gtk.Label(", "),
                        self.kpts_spin[2], gtk.Label(")")])
            self.kpts_label = gtk.Label("")
            self.kpts_label_format = _(u"k-points x size:  (%.1f, %.1f, %.1f) Å")
            pack(vbox, [self.kpts_label])
            self.k_changed()
            pack(vbox, gtk.Label(""))

        # Spin polarized, charge, relativity
        self.spinpol = gtk.CheckButton(_("Spin / initial moment "))
        self.spinpol.connect('toggled',self.spinpol_changed)
        self.moment  = gtk.Adjustment(0,-100,100,0.1)
        self.moment_spin = gtk.SpinButton(self.moment, 0, 0)
        self.moment_spin.set_digits(2)
        self.moment_spin.set_sensitive(False)
        self.charge  = gtk.Adjustment(0,-100,100,0.1)
        self.charge_spin = gtk.SpinButton(self.charge, 0, 0)
        self.charge_spin.set_digits(2)
        self.relativity_type = gtk.combo_box_new_text()
        for i, x in enumerate(self.aims_relativity_list):
            self.relativity_type.append_text(x)
        self.relativity_type.connect('changed',self.relativity_changed)
        self.relativity_threshold = gtk.Entry(max=8)
        self.relativity_threshold.set_text('1.00e-12')
        self.relativity_threshold.set_sensitive(False)
        pack(vbox, [self.spinpol,
                    self.moment_spin, 
                    gtk.Label(_("   Charge")), 
                    self.charge_spin, 
                    gtk.Label(_("   Relativity")),
                    self.relativity_type,
                    gtk.Label(_(" Threshold")),
                    self.relativity_threshold])
        pack(vbox, gtk.Label(""))

        # self-consistency criteria
        pack(vbox,[gtk.Label(_("Self-consistency convergence:"))])
        self.sc_tot_energy      = gtk.Adjustment(1e-6, 1e-6, 1e0, 1e-6)
        self.sc_tot_energy_spin = gtk.SpinButton(self.sc_tot_energy, 0, 0)
        self.sc_tot_energy_spin.set_digits(6)
        self.sc_tot_energy_spin.set_numeric(True)
        self.sc_sum_eigenvalue      = gtk.Adjustment(1e-3, 1e-6, 1e0, 1e-6)
        self.sc_sum_eigenvalue_spin = gtk.SpinButton(self.sc_sum_eigenvalue, 0, 0)
        self.sc_sum_eigenvalue_spin.set_digits(6)
        self.sc_sum_eigenvalue_spin.set_numeric(True)
        self.sc_density      = gtk.Adjustment(1e-4, 1e-6, 1e0, 1e-6)
        self.sc_density_spin = gtk.SpinButton(self.sc_density, 0, 0)
        self.sc_density_spin.set_digits(6)
        self.sc_density_spin.set_numeric(True)
        self.compute_forces = gtk.CheckButton(_("Compute forces"))
        self.compute_forces.set_active(True)
        self.compute_forces.connect("toggled", self.compute_forces_toggled,"")
        self.sc_forces      = gtk.Adjustment(1e-4, 1e-6, 1e0, 1e-6)
        self.sc_forces_spin = gtk.SpinButton(self.sc_forces, 0, 0)
        self.sc_forces_spin.set_numeric(True)
        self.sc_forces_spin.set_digits(6)
        # XXX: use gtk table for layout.  Spaces will not work well otherwise
        # (depend on fonts, widget style, ...)
        # TRANSLATORS: Don't care too much about these, just get approximately
        # the same string lengths
        pack(vbox, [gtk.Label(_("Energy:                 ")),
                    self.sc_tot_energy_spin, 
                    gtk.Label(_(" eV   Sum of eigenvalues:  ")),
                    self.sc_sum_eigenvalue_spin,
                    gtk.Label(_(" eV"))])
        pack(vbox, [gtk.Label(_("Electron density: ")),
                    self.sc_density_spin,
                    gtk.Label(_("        Force convergence:  ")),
                    self.sc_forces_spin,
                    gtk.Label(_(" eV/Ang  "))])

        pack(vbox, [self.compute_forces])
        pack(vbox, gtk.Label(""))

        swin = gtk.ScrolledWindow()
        swin.set_border_width(0)
        swin.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)

        self.expert_keyword_set = gtk.Entry(max = 55)
        self.expert_keyword_add = gtk.Button(stock = gtk.STOCK_ADD)
        self.expert_keyword_add.connect("clicked", self.expert_keyword_import)
        self.expert_keyword_set.connect("activate", self.expert_keyword_import)
        pack(vbox,[gtk.Label(_("Additional keywords: ")),
                   self.expert_keyword_set, 
                   self.expert_keyword_add])

        self.expert_vbox = gtk.VBox()
        vbox.pack_start(swin, True, True, 0)
        swin.add_with_viewport(self.expert_vbox)
        self.expert_vbox.get_parent().set_shadow_type(gtk.SHADOW_NONE)
        self.expert_vbox.get_parent().set_size_request(-1, 100)
        swin.show()
        self.expert_vbox.show()
        pack(vbox, gtk.Label(""))

        # run command and species defaults:
        pack(vbox, gtk.Label(_('FHI-aims execution command: ')))
        self.run_command = pack(vbox, gtk.Entry(max=0))
        pack(vbox, gtk.Label(_('Directory for species defaults: ')))
        self.species_defaults = pack(vbox, gtk.Entry(max=0))

        # set defaults from previous instance of the calculator, if applicable:
        if param is not None:
            self.set_param(param)
        else:
            self.set_defaults()

        # Buttons at the bottom
        pack(vbox, gtk.Label(""))
        butbox = gtk.HButtonBox()
        default_but = gtk.Button(_("Set Defaults"))
        default_but.connect("clicked",self.set_defaults)
        import_control_but = gtk.Button(_("Import control.in"))
        import_control_but.connect("clicked",self.import_control)
        export_control_but = gtk.Button(_("Export control.in"))
        export_control_but.connect("clicked", self.export_control)
        cancel_but = gtk.Button(stock=gtk.STOCK_CANCEL)
        cancel_but.connect('clicked', lambda widget: self.destroy())
        ok_but = gtk.Button(stock=gtk.STOCK_OK)
        ok_but.connect('clicked', self.ok)
        butbox.pack_start(default_but, 0, 0)
        butbox.pack_start(import_control_but, 0, 0)
        butbox.pack_start(export_control_but, 0, 0)
        butbox.pack_start(cancel_but, 0, 0)
        butbox.pack_start(ok_but, 0, 0)
        butbox.show_all()
        pack(vbox, [butbox], end=True, bottom=True)
        self.expert_vbox.show()
        vbox.show()
        self.add(vbox)
        self.show()
        self.grab_add() 
        if aims_periodic_warning:
            oops(aims_pbc_warning_text)

    def set_defaults(self, *args):
        atoms = self.owner.atoms.copy()
        if not self.xc_setup:
            self.xc_setup = True 
            for i, x in enumerate(self.xc_list):
                self.xc.append_text(x)
        for i, x in enumerate(self.xc_list):
            if x == self.aims_xc_default:
                self.xc.set_active(i)
        self.TS.set_active(False)
        if self.periodic:
            self.ucell = atoms.get_cell()
            for i in range(3):
                default = np.ceil(20.0 / np.sqrt(np.vdot(self.ucell[i],self.ucell[i])))
                self.kpts_spin[i].set_value(default)
        self.spinpol.set_active(False)
        self.moment.set_value(0)
        self.moment_spin.set_sensitive(False)
        self.charge.set_value(0)
        aims_relativity_default = 'none'
        for a in atoms:
            if a.number > 20: 
                aims_relativity_default = 'atomic_zora'
        for i, x in enumerate(self.aims_relativity_list):
            if x == aims_relativity_default:
                self.relativity_type.set_active(i)
        self.sc_tot_energy.set_value(1e-6)
        self.sc_sum_eigenvalue.set_value(1e-3)
        self.sc_density.set_value(1e-4)
        self.sc_forces.set_value(1e-4)
        for key in self.expert_keywords:
            key[0].destroy()
            key[1].destroy()
            key[2].destroy()
            key[3] = False
        for child in self.expert_vbox.children():
            self.expert_vbox.remove(child)
        if os.environ.has_key('AIMS_COMMAND'):
            text = os.environ['AIMS_COMMAND']
        else:
            text = ""
        self.run_command.set_text(text)
        if os.environ.has_key('AIMS_SPECIES_DIR'):
            text = os.environ['AIMS_SPECIES_DIR']
        else:
            text = ""
        self.species_defaults.set_text(text)

    def set_attributes(self, *args):
        param = {}
        param["xc"] = self.xc.get_active_text()
        if self.periodic:
            param["k_grid"] = (int(self.kpts[0].value),
                               int(self.kpts[1].value),
                               int(self.kpts[2].value))
        if self.spinpol.get_active():
            param["spin"] = "collinear"
            param["default_initial_moment"] = self.moment.get_value()
        else:
            param["spin"] = "none"
            param["default_initial_moment"] = None
        param["vdw_correction_hirshfeld"] = self.TS.get_active()
        param["charge"]             = self.charge.value
        param["relativistic"]       = self.relativity_type.get_active_text()
        if param["relativistic"] == 'atomic_zora':
            param["relativistic"] += " scalar "
        if param["relativistic"] == 'zora':
            param["relativistic"] += " scalar "+self.relativity_threshold.get_text() 
        param["sc_accuracy_etot"]   = self.sc_tot_energy.value
        param["sc_accuracy_eev"]    = self.sc_sum_eigenvalue.value
        param["sc_accuracy_rho"]    = self.sc_density.value
        param["compute_forces"]     = self.compute_forces.get_active()
        param["sc_accuracy_forces"] = self.sc_forces.value
        param["run_command"]        = self.run_command.get_text()
        param["species_dir"]        = self.species_defaults.get_text()
        from ase.calculators.aims import float_keys,exp_keys,string_keys,int_keys,bool_keys,list_keys,input_keys
        for option in self.expert_keywords:
            if option[3]:   # set type of parameter according to which list it is in
                key = option[0].get_text().strip()
                val = option[1].get_text().strip()
                if key == 'output':
                    if param.has_key('output'): 
                        param[key] += [val]
                    else:
                        param[key] = [val]
                elif key in float_keys or key in exp_keys:
                    param[key] = float(val)
                elif key in list_keys or key in string_keys or key in input_keys:
                    param[key] = val
                elif key in int_keys:
                    param[key] = int(val)
                elif key in bool_keys:
                    param[key] = bool(val)
        setattr(self.owner, self.attrname, param)

    def set_param(self, param):
        if param["xc"] is not None:
            for i, x in enumerate(self.xc_list):
                if x == param["xc"]:
                    self.xc.set_active(i)
        if isinstance(param["vdw_correction_hirshfeld"],bool):
            self.TS.set_active(param["vdw_correction_hirshfeld"])
        if self.periodic and param["k_grid"] is not None:
            self.kpts[0].value = int(param["k_grid"][0])
            self.kpts[1].value = int(param["k_grid"][1])
            self.kpts[2].value = int(param["k_grid"][2])
        if param["spin"] is not None:
            self.spinpol.set_active(param["spin"] == "collinear")
            self.moment_spin.set_sensitive(param["spin"] == "collinear")
        if param["default_initial_moment"] is not None:
            self.moment.value = param["default_initial_moment"]
        if param["charge"] is not None:
            self.charge.value = param["charge"]
        if param["relativistic"] is not None:
            if isinstance(param["relativistic"],(tuple,list)):
                rel = param["relativistic"]
            else:
                rel = param["relativistic"].split()
            for i, x in enumerate(self.aims_relativity_list):
                if x == rel[0]:
                    self.relativity_type.set_active(i)
                    if x == 'zora':
                        self.relativity_threshold.set_text(rel[2])
                        self.relativity_threshold.set_sensitive(True)
        if param["sc_accuracy_etot"] is not None:
            self.sc_tot_energy.value     = param["sc_accuracy_etot"]
        if param["sc_accuracy_eev"] is not None:
            self.sc_sum_eigenvalue.value = param["sc_accuracy_eev"]
        if param["sc_accuracy_rho"] is not None:
            self.sc_density.value        = param["sc_accuracy_rho"]
        if param["compute_forces"] is not None:
            if param["compute_forces"]:
                if param["sc_accuracy_forces"] is not None:
                    self.sc_forces.value = param["sc_accuracy_forces"]
                self.compute_forces.set_active(param["compute_forces"])
            else: 
                self.compute_forces.set_active(False)
        if param["run_command"] is not None:
            self.run_command.set_text(param["run_command"])
        if param["species_dir"] is not None:
            self.species_defaults.set_text(param["species_dir"])
        for (key,val) in param.items():
            if key in self.aims_keyword_list and key not in self.aims_keyword_gui_list:
                if val is not None:  # = existing "expert keyword"
                    if key == 'output': # 'output' can be used more than once
                        options = val
                        if isinstance(options,str): options = [options]
                        for arg in options:
                            self.expert_keyword_create([key]+[arg])
                    else:
                        if isinstance(val,str):
                            arg = [key]+val.split()
                        elif isinstance(val,(tuple,list)):
                            arg = [key]+[str(a) for a in val]
                        else:
                            arg = [key]+[str(val)]
                        self.expert_keyword_create(arg)

    def ok(self, *args):
        self.set_attributes(*args)
        self.destroy()

    def export_control(self, *args):
        filename = "control.in"
        chooser = gtk.FileChooserDialog(
            _('Export parameters ... '), None, gtk.FILE_CHOOSER_ACTION_SAVE,
            (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
             gtk.STOCK_SAVE, gtk.RESPONSE_OK))
        chooser.set_filename(filename)
        save = chooser.run()
        if save == gtk.RESPONSE_OK or save == gtk.RESPONSE_SAVE:
            filename = chooser.get_filename()
            self.set_attributes(*args)
            param = getattr(self.owner, "aims_parameters")
            from ase.calculators.aims import Aims
            calc_temp = Aims(**param)
            atoms_temp = self.owner.atoms.copy()
            atoms_temp.set_calculator(calc_temp)
            atoms_temp.calc.write_control(file = filename)
            atoms_temp.calc.write_species(file = filename)
        chooser.destroy()

    def import_control(self, *args):
        filename = "control.in"
        chooser = gtk.FileChooserDialog(
            _('Import control.in file ... '), None, 
            gtk.FILE_CHOOSER_ACTION_SAVE,
            (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
             gtk.STOCK_SAVE, gtk.RESPONSE_OK))
        chooser.set_filename(filename)
        save = chooser.run()
        if save == gtk.RESPONSE_OK:
            self.set_defaults()
            filename = chooser.get_filename()
            control = open(filename,'r')
            while True:
                line = control.readline()
                if not line:
                    break
                if "List of parameters used to initialize the calculator:" in line:
                    control.readline()
                    from ase.io.aims import read_aims_calculator
                    calc = read_aims_calculator(control)
                    found_aims_calculator = True
            control.close()
            if found_aims_calculator:
                param = calc.float_params
                for key in calc.exp_params:
                    param[key] = calc.exp_params[key]
                for key in calc.string_params:
                    param[key] = calc.string_params[key]
                for key in calc.int_params:
                    param[key] = calc.int_params[key]
                for key in calc.bool_params:
                    param[key] = calc.bool_params[key]
                for key in calc.list_params:
                    param[key] = calc.list_params[key]
                for key in calc.input_parameters:
                    param[key] = calc.input_parameters[key]
                self.set_defaults()
                self.set_param(param)
        chooser.destroy()

    def k_changed(self, *args):
        size = [self.kpts[i].value * np.sqrt(np.vdot(self.ucell[i],self.ucell[i])) for i in range(3)]
        self.kpts_label.set_text(self.kpts_label_format % tuple(size))

    def compute_forces_toggled(self, *args):
        self.sc_forces_spin.set_sensitive(self.compute_forces.get_active())

    def relativity_changed(self, *args):
        self.relativity_threshold.set_sensitive(self.relativity_type.get_active() == 2)

    def spinpol_changed(self, *args):
        self.moment_spin.set_sensitive(self.spinpol.get_active())

    def expert_keyword_import(self, *args):
        command = self.expert_keyword_set.get_text().split()
        if len(command) > 0 and command[0] in self.aims_keyword_list and not command[0] in self.aims_keyword_gui_list:
            self.expert_keyword_create(command)
        elif command[0] in self.aims_keyword_gui_list:
            oops(_("Please use the facilities provided in this window to "
                   "manipulate the keyword: %s!") % command[0])
        else:
            oops(_("Don't know this keyword: %s\n"
                   "\nPlease check!\n\n" 
                   "If you really think it should be available, "
                   "please add it to the top of ase/calculators/aims.py.")
                 % command[0])
        self.expert_keyword_set.set_text("")

    def expert_keyword_create(self, command):
        key = command[0]
        argument = command[1]
        if len(command) > 2:
            for a in command[2:]:
                argument += ' '+a
        index = len(self.expert_keywords) 
        self.expert_keywords += [[gtk.Label("    " +key+"  "),
                                  gtk.Entry(max=45),
                                  ExpertDeleteButton(index),
                                  True]]
        self.expert_keywords[index][1].set_text(argument)
        self.expert_keywords[index][2].connect('clicked',self.expert_keyword_delete)
        if not self.expert_vbox.get_children():
            table = gtk.Table(1, 3)
            table.attach(self.expert_keywords[index][0], 0, 1, 0, 1, 0)
            table.attach(self.expert_keywords[index][1], 1, 2, 0, 1, 0)
            table.attach(self.expert_keywords[index][2], 2, 3, 0, 1, 0)
            table.show_all()
            pack(self.expert_vbox, table)
        else:
            table = self.expert_vbox.get_children()[0]
            nrows = table.get_property('n-rows')
            table.resize(nrows + 1, 3)
            table.attach(self.expert_keywords[index][0],  0, 1, nrows, nrows + 1, 0) 
            table.attach(self.expert_keywords[index][1],  1, 2, nrows, nrows + 1, 0) 
            table.attach(self.expert_keywords[index][2],  2, 3, nrows, nrows + 1, 0) 
            table.show_all()

    def expert_keyword_delete(self, button, *args):
        index = button.index   # which one to kill 
        for i in [0,1,2]:
            self.expert_keywords[index][i].destroy()
        table = self.expert_vbox.get_children()[0]
        nrows = table.get_property('n-rows')
        table.resize(nrows-1, 3)
        self.expert_keywords[index][3] = False


class ExpertDeleteButton(gtk.Button):
    def __init__(self, index):
        gtk.Button.__init__(self, stock=gtk.STOCK_DELETE)
        alignment = self.get_children()[0]
        hbox = alignment.get_children()[0]
        #self.set_size_request(1, 3)
        image, label = hbox.get_children()
        if image is not None: 
            label.set_text('Del')
        self.index = index


class VASP_Window(gtk.Window):
    vasp_xc_list = ['PW91', 'PBE', 'LDA']
    vasp_xc_default = 'PBE'
    vasp_prec_default = 'Normal'
    def __init__(self, owner, param, attrname):
        self.owner = owner
        self.attrname = attrname
        atoms = owner.atoms
        self.periodic = atoms.get_pbc().all()
        self.vasp_keyword_gui_list = ['ediff','encut', 'ismear', 'ispin', 'prec', 'sigma']
        from ase.calculators.vasp import float_keys,exp_keys,string_keys,int_keys,bool_keys,list_keys,special_keys
        self.vasp_keyword_list = float_keys+exp_keys+string_keys+int_keys+bool_keys+list_keys+special_keys
        self.expert_keywords = []
        natoms = len(atoms)
        gtk.Window.__init__(self)
        self.set_title(_("VASP parameters"))
        vbox = gtk.VBox()
        vbox.set_border_width(5)
        # Print some info
        txt = _("%i atoms.\n") % natoms
        self.ucell = atoms.get_cell()
        txt += _("Periodic geometry, unit cell is: \n")
        for i in range(3):
            txt += "(%8.3f %8.3f %8.3f)\n" % (self.ucell[i][0], self.ucell[i][1], self.ucell[i][2])
        pack(vbox, [gtk.Label(txt)])

        # XC functional ()
        self.xc = gtk.combo_box_new_text()
        for i, x in enumerate(self.vasp_xc_list):
            self.xc.append_text(x)

        # Spin polarized
        self.spinpol = gtk.CheckButton(_("Spin polarized"))
        
        pack(vbox, [gtk.Label(_("Exchange-correlation functional: ")),
                    self.xc,
                    gtk.Label("    "),
                    self.spinpol])
        pack(vbox, gtk.Label(""))

        # k-grid
        self.kpts = []
        self.kpts_spin = []
        for i in range(3):
            default = np.ceil(20.0 / np.sqrt(np.vdot(self.ucell[i],self.ucell[i])))
            g = gtk.Adjustment(default, 1, 100, 1)
            s = gtk.SpinButton(g, 0, 0)
            self.kpts.append(g)
            self.kpts_spin.append(s)
            g.connect("value-changed", self.k_changed)

        # Precision of calculation
        self.prec = gtk.combo_box_new_text()
        for i, x in enumerate(['Low', 'Normal', 'Accurate']):
            self.prec.append_text(x)
            if x == self.vasp_prec_default:
                self.prec.set_active(i)

        # cutoff energy
        if os.environ.has_key('VASP_PP_PATH'):
            self.encut_min_default, self.encut_max_default = self.get_min_max_cutoff()
        else:
            self.encut_max_default = 400.0
            self.encut_min_default = 100.0
        self.encut = gtk.Adjustment(self.encut_max_default, 0, 9999, 10)
        self.encut_spin = gtk.SpinButton(self.encut, 0, 0)
        self.encut_spin.set_digits(2)
        self.encut_spin.connect("value-changed",self.check_encut_warning)
        self.encut_warning = gtk.Label("")

        pack(vbox, [gtk.Label(_("k-points  k = (")), self.kpts_spin[0],
                    gtk.Label(", "), self.kpts_spin[1], gtk.Label(", "),
                    self.kpts_spin[2], 
                    gtk.Label(_(")    Cutoff: ")),self.encut_spin,
                    gtk.Label(_("    Precision: ")),self.prec])
        self.kpts_label = gtk.Label("")
        self.kpts_label_format = _(u"k-points x size:  (%.1f, %.1f, %.1f) Å       ")
        pack(vbox, [self.kpts_label, self.encut_warning])
        self.k_changed()
        pack(vbox, gtk.Label(""))

        self.ismear = gtk.combo_box_new_text()
        for x in ['Fermi', 'Gauss', 'Methfessel-Paxton']:
            self.ismear.append_text(x)
        self.ismear.set_active(2)
        self.smearing_order = gtk.Adjustment(2,0,9,1)
        self.smearing_order_spin = gtk.SpinButton(self.smearing_order,0,0)
        self.smearing_order_spin.set_digits(0)
        self.ismear.connect("changed", self.check_ismear_changed)
        self.sigma = gtk.Adjustment(0.1, 0.001, 9.0, 0.1)
        self.sigma_spin = gtk.SpinButton(self.sigma,0,0)
        self.sigma_spin.set_digits(3)
        pack(vbox, [gtk.Label(_("Smearing: ")),
                    self.ismear,
                    gtk.Label(_(" order: ")),
                    self.smearing_order_spin,
                    gtk.Label(_(" width: ")),
                    self.sigma_spin])
        pack(vbox, gtk.Label(""))
        
        self.ediff = gtk.Adjustment(1e-4, 1e-6, 1e0, 1e-4)
        self.ediff_spin = gtk.SpinButton(self.ediff, 0, 0)
        self.ediff_spin.set_digits(6)
        pack(vbox,[gtk.Label(_("Self-consistency convergence: ")),
                   self.ediff_spin,
                   gtk.Label(_(" eV"))])
        pack(vbox,gtk.Label(""))

        swin = gtk.ScrolledWindow()
        swin.set_border_width(0)
        swin.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)

        self.expert_keyword_set = gtk.Entry(max = 55)
        self.expert_keyword_add = gtk.Button(stock = gtk.STOCK_ADD)
        self.expert_keyword_add.connect("clicked", self.expert_keyword_import)
        self.expert_keyword_set.connect("activate", self.expert_keyword_import)
        pack(vbox,[gtk.Label(_("Additional keywords: ")),
                   self.expert_keyword_set, 
                   self.expert_keyword_add])
        self.expert_vbox = gtk.VBox()
        vbox.pack_start(swin, True, True, 0)
        swin.add_with_viewport(self.expert_vbox)
        self.expert_vbox.get_parent().set_shadow_type(gtk.SHADOW_NONE)
        self.expert_vbox.get_parent().set_size_request(-1, 100)
        swin.show()
        self.expert_vbox.show()
        pack(vbox, gtk.Label(""))

        # run command and location of POTCAR files:
        pack(vbox, gtk.Label(_('VASP execution command: ')))
        self.run_command = pack(vbox, gtk.Entry(max=0))
        if os.environ.has_key('VASP_COMMAND'):
            self.run_command.set_text(os.environ['VASP_COMMAND'])
        pack(vbox, gtk.Label(_('Directory for species defaults: ')))
        self.pp_path = pack(vbox, gtk.Entry(max=0))
        if os.environ.has_key('VASP_PP_PATH'):
            self.pp_path.set_text(os.environ['VASP_PP_PATH'])

        # Buttons at the bottom
        pack(vbox, gtk.Label(""))
        butbox = gtk.HButtonBox()
        set_default_but = gtk.Button(_("Set Defaults"))
        set_default_but.connect("clicked", self.set_defaults)
        import_vasp_but = gtk.Button(_("Import VASP files"))
        import_vasp_but.connect("clicked", self.import_vasp_files)
        export_vasp_but = gtk.Button(_("Export VASP files"))
        export_vasp_but.connect("clicked", self.export_vasp_files)
        cancel_but = gtk.Button(stock=gtk.STOCK_CANCEL)
        cancel_but.connect('clicked', lambda widget: self.destroy())
        ok_but = gtk.Button(stock=gtk.STOCK_OK)
        ok_but.connect('clicked', self.ok)
        butbox.pack_start(set_default_but, 0, 0)
        butbox.pack_start(import_vasp_but, 0, 0)
        butbox.pack_start(export_vasp_but, 0, 0)
        butbox.pack_start(cancel_but, 0, 0)
        butbox.pack_start(ok_but, 0, 0)
        butbox.show_all()
        pack(vbox, [butbox], end=True, bottom=True)
        vbox.show()
        self.add(vbox)
        self.show()
        self.grab_add()  # Lock all other windows

        self.load_attributes()

    def load_attributes(self, directory = "."):
        """Sets values of fields of the window according to the values 
        set inside the INCAR, KPOINTS and POTCAR file in 'directory'."""
        from os import chdir
        chdir(directory)
       
        # Try and load INCAR, in the current directory
        from ase.calculators.vasp import Vasp
        calc_temp = Vasp()
        try:
            calc_temp.read_incar("INCAR")
        except IOError:
            pass
        else:
            if calc_temp.spinpol:
                self.spinpol.set_active(True)
            else:
                self.spinpol.set_active(False)

            if calc_temp.float_params['encut']:
                self.encut.set_value(calc_temp.float_params['encut'])
 
            if calc_temp.int_params['ismear'] == -1: # Fermi
                vasp_ismear_default = 'Fermi'
            elif calc_temp.int_params['ismear'] == 0: # Gauss
                vasp_ismear_default = 'Gauss'
            elif calc_temp.int_params['ismear'] > 0: # Methfessel-Paxton
                vasp_ismear_default = 'Methfessel-Paxton'
            else:
                vasp_ismear_default = None

            for i, x in enumerate(['Fermi', 'Gauss', 'Methfessel-Paxton']):
                if vasp_ismear_default == x:
                    self.ismear.set_active(i)

            if calc_temp.exp_params['ediff']:
                self.ediff.set_value(calc_temp.exp_params['ediff'])

            for i, x in enumerate(['Low', 'Normal', 'Accurate']):
                if x == calc_temp.string_params['prec']:
                    self.prec.set_active(i)

            if calc_temp.float_params['sigma']:
                self.sigma.set_value(calc_temp.float_params['sigma'])

            import copy
            all_params = copy.deepcopy(calc_temp.float_params)
            all_params.update(calc_temp.exp_params)
            all_params.update(calc_temp.string_params)
            all_params.update(calc_temp.int_params)
            all_params.update(calc_temp.bool_params)
            all_params.update(calc_temp.special_params)

            for (key, value) in all_params.items(): 
                if key in self.vasp_keyword_list \
                        and key not in self.vasp_keyword_gui_list \
                        and value is not None:
                    command = key + " " + str(value)
                    self.expert_keyword_create(command.split())

            for (key, value) in calc_temp.list_params.items():
                if key == "magmom" and value is not None:
                    command = key + " "
                    rep = 1
                    previous = value[0]
                    for v in value[1:]:
                        if v == previous:
                            rep += 1
                        else:
                            if rep > 1:
                                command += "%d*%f " % (rep, previous)
                            else:
                                command += "%f " % previous
                            rep = 1
                        previous = v
                    if rep > 1:
                        command += "%d*%f " % (rep, previous)
                    else:
                        command += "%f" % previous
                    self.expert_keyword_create(command.split())
                elif value is not None:
                    command = key + " "
                    for v in value:
                        command += str(v) + " "
                    self.expert_keyword_create(command.split())
                 


        # Try and load POTCAR, in the current directory
        try:
            calc_temp.read_potcar()
        except IOError:
            pass
        else:
            #Set xc read from POTCAR
            for i, x in enumerate(self.vasp_xc_list):
                if x == calc_temp.input_params['xc']:
                    self.xc.set_active(i)

        # Try and load KPOINTS, in the current directory
        try:
            calc_temp.read_kpoints("KPOINTS")
        except IOError:
            pass
        else:
            # Set KPOINTS grid dimensions
            for i in range(3):
                self.kpts_spin[i].set_value(calc_temp.input_params['kpts'][i])

    def set_attributes(self, *args):
        self.param = {}
        self.param["xc"] = self.xc.get_active_text()
        self.param["prec"] = self.prec.get_active_text()
        self.param["kpts"] = (int(self.kpts[0].value),
                              int(self.kpts[1].value),
                              int(self.kpts[2].value))
        self.param["encut"] = self.encut.value
        self.param["ediff"] = self.ediff.value
        self.param["ismear"] = self.get_ismear()
        self.param["sigma"] = self.sigma.value
        if self.spinpol.get_active(): 
            self.param["ispin"] = 2
        else:
            self.param["ispin"] = 1
        from ase.calculators.vasp import float_keys,exp_keys,string_keys,int_keys,bool_keys,list_keys,special_keys
        for option in self.expert_keywords:
            if option[3]:   # set type of parameter accoding to which list it is in
                key = option[0].get_text().split()[0].strip()
                val = option[1].get_text().strip()
                if key in float_keys or key in exp_keys:
                    self.param[key] = float(val)
                elif key == "magmom":
                    val = val.replace("*", " * ")
                    c = val.split()
                    val = []
                    i = 0
                    while i < len(c):
                        if c[i] == "*":
                            b = val.pop()
                            i += 1
                            for j in range(int(b)):
                                val.append(float(c[i]))
                        else:
                            val.append(float(c[i]))
                        i += 1
                    self.param[key] = val
                elif key in list_keys:
                    c = val.split()
                    val = []
                    for i in c:
                        val.append(float(i))
                    self.param[key] = val
                elif key in string_keys or key in special_keys:
                    self.param[key] = val
                elif key in int_keys:
                    self.param[key] = int(val)
                elif key in bool_keys:
                    self.param[key] = bool(val)
        setattr(self.owner, self.attrname, self.param)
        os.environ['VASP_COMMAND'] = self.run_command.get_text()
        os.environ['VASP_PP_PATH'] = self.pp_path.get_text()
        
    def ok(self, *args):
        self.set_attributes(*args)
        self.destroy()

    def get_min_max_cutoff(self, *args):
        # determine the recommended energy cutoff limits 
        from ase.calculators.vasp import Vasp
        calc_temp = Vasp()
        atoms_temp = self.owner.atoms.copy()
        calc_temp.initialize(atoms_temp)
        calc_temp.write_potcar(suffix = '.check_energy_cutoff')
        enmin = -1e6
        enmax = -1e6
        for line in open("POTCAR.check_energy_cutoff",'r').readlines():
            if "ENMIN" in line:
                enmax = max(enmax,float(line.split()[2].split(';')[0]))
                enmin = max(enmin,float(line.split()[5]))
        from os import system
        system("rm POTCAR.check_energy_cutoff")
        return enmin, enmax

    def k_changed(self, *args):
        size = [self.kpts[i].value * np.sqrt(np.vdot(self.ucell[i],self.ucell[i])) for i in range(3)]
        self.kpts_label.set_text(self.kpts_label_format % tuple(size))

    def check_encut_warning(self,*args):
        if self.encut.value < self.encut_min_default:
            self.encut_warning.set_markup(_("<b>WARNING:</b> cutoff energy is lower than recommended minimum!"))
        else:
            self.encut_warning.set_markup("")

    def check_ismear_changed(self,*args):
        if self.ismear.get_active_text() == 'Methfessel-Paxton':
            self.smearing_order_spin.set_sensitive(True)
        else:
            self.smearing_order_spin.set_sensitive(False)

    def get_ismear(self,*args):
        type = self.ismear.get_active_text()
        if type == 'Methfessel-Paxton':
            ismear_value = self.smearing_order.value
        elif type == 'Fermi':
            ismear_value = -1
        else:
            ismear_value = 0
        return ismear_value

    def destroy(self):
        self.grab_remove()
        gtk.Window.destroy(self)

    def set_defaults(self, *args):
         # Reset fields to what they were
        self.spinpol.set_active(False)

        for i, x in enumerate(['Low', 'Normal', 'Accurate']):
            if x == self.vasp_prec_default:
                self.prec.set_active(i)

        self.encut_spin.set_value(self.encut_max_default)

        self.ismear.set_active(2)
        self.smearing_order.set_value(2)
        self.ediff.set_value(1e-4)

        for child in self.expert_vbox.children():
            self.expert_vbox.remove(child)

        for i, x in enumerate(self.vasp_xc_list):
                if x == self.vasp_xc_default:
                    self.xc.set_active(i) 

        default = np.ceil(20.0 / np.sqrt(np.vdot(self.ucell[i],self.ucell[i])))
        for i in range(3):
            self.kpts_spin[i].set_value(default)

    def import_vasp_files(self, *args):
        dirname = ""
        chooser = gtk.FileChooserDialog(
            _('Import VASP input files: choose directory ... '),
            None, gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER,
            (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
             gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        chooser.set_filename(dirname)
        openr = chooser.run()
        if openr == gtk.RESPONSE_OK or openr == gtk.RESPONSE_SAVE:
            dirname = chooser.get_filename()
            self.load_attributes(dirname)
        chooser.destroy()
            

    def export_vasp_files(self, *args):
        filename = ""
        chooser = gtk.FileChooserDialog(
            _('Export VASP input files: choose directory ... '), 
            None, gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER,
            (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
             gtk.STOCK_SAVE, gtk.RESPONSE_OK))
        chooser.set_filename(filename)
        save = chooser.run()
        if save == gtk.RESPONSE_OK or save == gtk.RESPONSE_SAVE:
            filename = chooser.get_filename()
            from os import chdir
            chdir(filename)
            self.set_attributes(*args)
            param = getattr(self.owner, "vasp_parameters")
            from ase.calculators.vasp import Vasp
            calc_temp = Vasp(**param)
            atoms_temp = self.owner.atoms.copy()
            atoms_temp.set_calculator(calc_temp)
            calc_temp.initialize(atoms_temp)
            calc_temp.write_incar(atoms_temp)
            calc_temp.write_potcar()
            calc_temp.write_kpoints()
            calc_temp.write_sort_file()
            from ase.io.vasp import write_vasp
            write_vasp('POSCAR', calc_temp.atoms_sorted, symbol_count = calc_temp.symbol_count)
        chooser.destroy()

    def expert_keyword_import(self, *args):
        command = self.expert_keyword_set.get_text().split()
        if len(command) > 0 and command[0] in self.vasp_keyword_list and not command[0] in self.vasp_keyword_gui_list:
            self.expert_keyword_create(command)
        elif command[0] in self.vasp_keyword_gui_list:
            oops(_("Please use the facilities provided in this window to "
                   "manipulate the keyword: %s!") % command[0])
        else:
            oops(_("Don't know this keyword: %s"
                   "\nPlease check!\n\n" 
                   "If you really think it should be available, "
                   "please add it to the top of ase/calculators/vasp.py.")
                 % command[0])
        self.expert_keyword_set.set_text("")

    def expert_keyword_create(self, command):
        key = command[0]
        if command[1] == "=":
            command.remove("=")
        argument = command[1]
        if len(command) > 2:
            for a in command[2:]:
                argument += ' '+a
        index = len(self.expert_keywords) 
        self.expert_keywords += [[gtk.Label("    " +key+" = "),
                                  gtk.Entry(max=55),
                                  ExpertDeleteButton(index),
                                  True]]
        self.expert_keywords[index][1].set_text(argument)
        self.expert_keywords[index][2].connect('clicked',self.expert_keyword_delete)
        if not self.expert_vbox.get_children():
            table = gtk.Table(1, 3)
            table.attach(self.expert_keywords[index][0], 0, 1, 0, 1, 0)
            table.attach(self.expert_keywords[index][1], 1, 2, 0, 1, 0)
            table.attach(self.expert_keywords[index][2], 2, 3, 0, 1, 0)
            table.show_all()
            pack(self.expert_vbox, table)
        else:
            table = self.expert_vbox.get_children()[0]
            nrows = table.get_property('n-rows')
            table.resize(nrows + 1, 3)
            table.attach(self.expert_keywords[index][0],  0, 1, nrows, nrows + 1, 0) 
            table.attach(self.expert_keywords[index][1],  1, 2, nrows, nrows + 1, 0) 
            table.attach(self.expert_keywords[index][2],  2, 3, nrows, nrows + 1, 0) 
            table.show_all()
        
    def expert_keyword_delete(self, button, *args):
        index = button.index   # which one to kill 
        for i in [0,1,2]:
            self.expert_keywords[index][i].destroy()
        table = self.expert_vbox.get_children()[0]
        nrows = table.get_property('n-rows')
        table.resize(nrows-1, 3)
        self.expert_keywords[index][3] = False
