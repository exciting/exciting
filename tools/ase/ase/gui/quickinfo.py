# encoding: utf-8

"Module for displaying information about the system."

import gtk
from gettext import gettext as _
from ase.gui.widgets import pack

singleimage = _("Single image loaded.")
multiimage = _("Image %d loaded (0 - %d).")
ucconst = _("Unit cell is fixed.")
ucvaries = _("Unit cell varies.")

format = _("""\
%s

Number of atoms: %d.

Unit cell:
  %8.3f  %8.3f  %8.3f
  %8.3f  %8.3f  %8.3f
  %8.3f  %8.3f  %8.3f
%s
""")

class QuickInfo(gtk.Window):
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.set_title(_("Quick Info"))
        vbox = gtk.VBox()
        images = gui.images
        if images.natoms < 1:
            txt = _("No atoms loaded.")
        else:
            (nimg, natoms, three) = images.P.shape
            assert three == 3
            img = gui.frame
            uc = images.A[img]
            if nimg > 1:
                equal = True
                for i in range(nimg):
                    equal = equal and (uc == images.A[i]).all()
                if equal:
                    uctxt = ucconst
                else:
                    uctxt = ucvaries
            else:
                uctxt = ""
            if nimg == 1:
                imgtxt = singleimage
            else:
                imgtxt = multiimage % (img, nimg-1)
            txt = format % ((imgtxt, natoms) + tuple(uc.flat) + (uctxt,))
        label = gtk.Label(txt)
        pack(vbox, [label])
        but = gtk.Button(stock=gtk.STOCK_CLOSE)
        but.connect('clicked', self.close)
        pack(vbox, [but], end=True)
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui

    def close(self, *args):
        self.destroy()

    
