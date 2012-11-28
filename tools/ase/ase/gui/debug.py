import sys

import gtk
from gettext import gettext as _
import numpy as np

class Debug(gtk.Window):
    def __init__(self, gui):
        self.gui = gui
        gtk.Window.__init__(self)
        self.set_title(_('Debug'))
        entry = gtk.Entry(200)
        self.add(entry)
        entry.connect('activate', self.enter, entry)
        entry.show()
        self.show()

    def enter(self, widget, entry):
        g = self.gui
        print >> sys.stderr, eval(entry.get_text())
    
