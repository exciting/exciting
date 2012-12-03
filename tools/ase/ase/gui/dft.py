#!/usr/bin/env python
from math import sqrt

import gtk
from gettext import gettext as _

from ase.gui.widgets import pack, Help


class DFT(gtk.Window):
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.set_title(_('DFT'))
        vbox = gtk.VBox()
        combo = gtk.combo_box_new_text()
        self.xcfuncs = 'None LDA PBE revPBE RPBE PW91 EXX PBE0'.split()
        for xc in self.xcfuncs:
            combo.append_text(xc)
        pack(vbox, [gtk.Label(_('XC-functional: ')), combo])

        button=radio(None,monkhorstpack)
        button=radio(button, special)
        pack(vbox, gtk.Label(_('Repeat atoms:')))
        self.kpts = [gtk.Adjustment(r, 1, 99, 1) for r in gui.atoms.repeat]
        pack(vbox, [gtk.SpinButton(r, 0, 0) for r in self.repeat])
        for r in self.repeat:
            r.connect('value-changed', self.change)

        close = pack(vbox, gtk.Button(_('Close')))
        close.connect('clicked', lambda widget: self.destroy())
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui

        xc = gui.atoms.dft.get('xc', 'None')
        combo.set_active(self.xcfuncs.index(xc))

    def selected(self, button):
        self.gui.atoms.dynamic = ~self.gui.atoms.selected
        self.gui.draw()

    def immobile(self, button):
        self.gui.atoms.set_dynamic()
        self.gui.draw()

    def clear(self, button):
        self.gui.atoms.dynamic[:] = True
        self.gui.draw()

