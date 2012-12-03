#!/usr/bin/env python
import gtk
from math import sqrt
from gettext import gettext as _

import numpy as np

from ase.gui.widgets import pack, Help


class Repeat(gtk.Window):
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.set_title(_('Repeat'))
        vbox = gtk.VBox()
        pack(vbox, gtk.Label(_('Repeat atoms:')))
        self.repeat = [gtk.Adjustment(r, 1, 9, 1) for r in gui.images.repeat]
        pack(vbox, [gtk.SpinButton(r, 0, 0) for r in self.repeat])
        for r in self.repeat:
            r.connect('value-changed', self.change)
        button = pack(vbox, gtk.Button(_('Set unit cell')))
        button.connect('clicked', self.set_unit_cell)
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui

    def change(self, adjustment):
        self.gui.images.repeat_images([int(r.value) for r in self.repeat])
        self.gui.repeat_colors([int(r.value) for r in self.repeat])
        self.gui.set_coordinates()
        return True
        
    def set_unit_cell(self, button):
        self.gui.images.A *= self.gui.images.repeat.reshape((3, 1))
        self.gui.images.E *= self.gui.images.repeat.prod()
        self.gui.images.repeat = np.ones(3, int)
        for r in self.repeat:
            r.value = 1
        self.gui.set_coordinates()
        
