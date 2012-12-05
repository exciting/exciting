#!/usr/bin/env python
from math import sqrt

import gtk
from gettext import gettext as _

from ase.gui.widgets import pack, Help


class Constraints(gtk.Window):
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.set_title(_('Constraints'))
        vbox = gtk.VBox()
        b = pack(vbox, [gtk.Button(_('Constrain')),
                        gtk.Label(_(' selected atoms'))])[0]
        b.connect('clicked', self.selected)
        b = pack(vbox, [gtk.Button(_('Constrain')),
                        gtk.Label(_(' immobile atoms:'))])[0]
        b.connect('clicked', self.immobile)
        b = pack(vbox, [gtk.Button(_('Unconstrain')),
                        gtk.Label(_(' selected atoms:'))])[0]
        b.connect('clicked', self.unconstrain)
        b = pack(vbox, gtk.Button(_('Clear constraints')))
        b.connect('clicked', self.clear)
        close = pack(vbox, gtk.Button(_('Close')))
        close.connect('clicked', lambda widget: self.destroy())
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui

    def selected(self, button):
        self.gui.images.dynamic[self.gui.images.selected] = False
        self.gui.draw()

    def unconstrain(self, button):
        self.gui.images.dynamic[self.gui.images.selected] = True
        self.gui.draw()
        
    def immobile(self, button):
        self.gui.images.set_dynamic()
        self.gui.draw()

    def clear(self, button):
        self.gui.images.dynamic[:] = True
        self.gui.draw()

