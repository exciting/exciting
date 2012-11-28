# encoding: utf-8

"Module for performing energy minimization."

import gtk
from gettext import gettext as _
from ase.gui.simulation import Simulation
from ase.gui.widgets import oops, pack, AseGuiCancelException
import ase
import ase.optimize
import numpy as np

class MinimizeMixin:
    minimizers = ('BFGS', 'BFGSLineSearch', 'LBFGS', 'LBFGSLineSearch', 'MDMin', 'FIRE')
    def make_minimize_gui(self, box):
        self.algo = gtk.combo_box_new_text()
        for m in self.minimizers:
            self.algo.append_text(m)
        self.algo.set_active(0)
        self.algo.connect('changed', self.min_algo_specific)
        pack(box, [gtk.Label(_("Algorithm: ")), self.algo])
        
        self.fmax = gtk.Adjustment(0.05, 0.00, 10.0, 0.01)
        self.fmax_spin = gtk.SpinButton(self.fmax, 0, 3)
        lbl = gtk.Label()
        lbl.set_markup(_("Convergence criterion: F<sub>max</sub> = "))
        pack(box, [lbl, self.fmax_spin])

        self.steps = gtk.Adjustment(100, 1, 1000000, 1)
        self.steps_spin = gtk.SpinButton(self.steps, 0, 0)
        pack(box, [gtk.Label(_("Max. number of steps: ")), self.steps_spin])

        # Special stuff for MDMin
        lbl = gtk.Label(_("Pseudo time step: "))
        self.mdmin_dt = gtk.Adjustment(0.05, 0.0, 10.0, 0.01)
        spin = gtk.SpinButton(self.mdmin_dt, 0, 3)
        self.mdmin_widgets = [lbl, spin]
        pack(box, self.mdmin_widgets)
        self.min_algo_specific()
        
    def min_algo_specific(self, *args):
        "SHow or hide algorithm-specific widgets."
        minimizer = self.minimizers[self.algo.get_active()]
        for w in self.mdmin_widgets:
            if minimizer == 'MDMin':
                w.show()
            else:
                w.hide()
        
class Minimize(Simulation, MinimizeMixin):
    "Window for performing energy minimization."
    
    def __init__(self, gui):
        Simulation.__init__(self, gui)
        self.set_title(_("Energy minimization"))
        
        vbox = gtk.VBox()
        self.packtext(vbox,
                      _("Minimize the energy with respect to the positions."))
        self.packimageselection(vbox)
        pack(vbox, gtk.Label(""))

        self.make_minimize_gui(vbox)
        
        pack(vbox, gtk.Label(""))
        self.status_label = gtk.Label("")
        pack(vbox, [self.status_label])
        self.makebutbox(vbox)
        vbox.show()
        self.add(vbox)
        self.show()
        self.gui.register_vulnerable(self)

    def run(self, *args):
        "User has pressed [Run]: run the minimization."
        if not self.setup_atoms():
            return
        fmax = self.fmax.value
        steps = self.steps.value
        mininame = self.minimizers[self.algo.get_active()]
        self.begin(mode="min", algo=mininame, fmax=fmax, steps=steps)
        algo = getattr(ase.optimize, mininame)
        try:
            logger_func = self.gui.simulation['progress'].get_logger_stream
        except (KeyError, AttributeError):
            logger = None
        else:
            logger = logger_func()  # Don't catch errors in the function.

        # Display status message
        self.status_label.set_text(_("Running ..."))
        self.status_label.modify_fg(gtk.STATE_NORMAL,
                                    gtk.gdk.color_parse('#AA0000'))
        while gtk.events_pending():
            gtk.main_iteration()

        self.prepare_store_atoms()
        if mininame == "MDMin":
            minimizer = algo(self.atoms, logfile=logger,
                             dt=self.mdmin_dt.value)
        else:
            minimizer = algo(self.atoms, logfile=logger)
        minimizer.attach(self.store_atoms)
        try:
            minimizer.run(fmax=fmax, steps=steps)
        except AseGuiCancelException:
            # Update display to reflect cancellation of simulation.
            self.status_label.set_text(_("Minimization CANCELLED after "
                                         "%i steps.")
                                       % (self.count_steps,))
            self.status_label.modify_fg(gtk.STATE_NORMAL,
                                        gtk.gdk.color_parse('#AA4000'))
        except MemoryError:
            self.status_label.set_text(_("Out of memory, consider using "
                                         "LBFGS instead"))
            self.status_label.modify_fg(gtk.STATE_NORMAL,
                                        gtk.gdk.color_parse('#AA4000'))
            
        else:
            # Update display to reflect succesful end of simulation.
            self.status_label.set_text(_("Minimization completed in %i steps.")
                                       % (self.count_steps,))
            self.status_label.modify_fg(gtk.STATE_NORMAL,
                                        gtk.gdk.color_parse('#007700'))
            
        self.end()
        if self.count_steps:
            # Notify other windows that atoms have changed.
            # This also notifies this window!
            self.gui.notify_vulnerable()

        # Open movie window and energy graph
        if self.gui.images.nimages > 1:
            self.gui.movie()
            assert not np.isnan(self.gui.images.E[0])
            if not self.gui.plot_graphs_newatoms():
                expr = 'i, e - E[-1]'            
                self.gui.plot_graphs(expr=expr)

    def notify_atoms_changed(self):
        "When atoms have changed, check for the number of images."
        self.setupimageselection()
        
