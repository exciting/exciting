#!/usr/bin/env python
import gtk
from ase.gui.widgets import pack
from gettext import gettext as _

class Settings(gtk.Window):
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.set_title('Settings')
        self.gui = gui
        vbox = gtk.VBox()

        # Constraints
        a = pack(vbox, gtk.Label())
        a.set_markup('<span size="larger" underline="single">'
                     '%s</span>' % _('Constraints:'))
        a, b = pack(vbox, [gtk.Button(_('Constrain')),
                           gtk.Label('/'),
                           gtk.Button(_('release')),                        
                           gtk.Label(_(' selected atoms'))])[::2]
        a.connect('clicked', self.constrain_selected)
        b.connect('clicked', self.release_selected)
        a = pack(vbox, gtk.Button(_('Constrain immobile atoms')))
        a.connect('clicked', self.immobile)
        a = pack(vbox, gtk.Button(_('Clear all constraints')))
        a.connect('clicked', self.clear_constraints)

        # Visibility
        a = pack(vbox, gtk.Label())
        a.set_markup('\n<span size="larger" underline="single">'
                     '%s</span>' % _('Visibility:'))
        a, b = pack(vbox, [gtk.Button(_('Hide')),
                           gtk.Label('/'),
                           gtk.Button(_('show')),
                           gtk.Label(_(' selected atoms'))])[::2]
        a.connect('clicked', self.hide_selected)
        b.connect('clicked', self.show_selected)
        a = pack(vbox, gtk.Button(_('View all atoms')))
        a.connect('clicked', self.view_all)

        # Miscellaneous
        a = pack(vbox, gtk.Label())
        a.set_markup('\n<span size="larger" underline="single">'
                     '%s</span>' % _('Miscellaneous:'))
        self.scale = gtk.Adjustment(value=.89, lower=0.2, upper=2.0,
                                    step_incr=0.1, page_incr=0.5)
        pack(vbox, [gtk.Label(_('Scale atomic radii:')),
                    gtk.SpinButton(self.scale, climb_rate=0, digits=2)])
        self.scale.connect('value-changed', self.scale_radii)

        # A close button
        pack(vbox, gtk.Label(_('\n')))
        close = pack(vbox, gtk.Button(_('Close')))
        close.connect('clicked', lambda widget: self.destroy())

        # Add elements and show frame
        self.add(vbox)
        vbox.show()
        self.show()

    def scale_radii(self, adjustment):
        self.gui.images.set_radii(float(self.scale.value))
        self.gui.draw()
        return True

    def hide_selected(self, button):
        self.gui.images.visible[self.gui.images.selected] = False
        self.gui.draw()

    def show_selected(self, button):
        self.gui.images.visible[self.gui.images.selected] = True
        self.gui.draw()

    def view_all(self, button):
        self.gui.images.visible[:] = True
        self.gui.draw()

    def constrain_selected(self, button):
        self.gui.images.dynamic[self.gui.images.selected] = False
        self.gui.draw()

    def release_selected(self, button):
        self.gui.images.dynamic[self.gui.images.selected] = True
        self.gui.draw()

    def immobile(self, button):
        self.gui.images.set_dynamic()
        self.gui.draw()

    def clear_constraints(self, button):
        self.gui.images.dynamic[:] = True
        self.gui.draw()
