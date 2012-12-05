# -*- coding: utf-8 -*-

from math import sqrt, pi, acos

import gtk
import numpy as np

from ase.data import chemical_symbols as symbols
from ase.data import atomic_names as names
from ase.gui.widgets import pack
from gettext import gettext as _

def formula(Z):
    hist = {}
    for z in Z:
        if z in hist:
            hist[z] += 1
        else:
            hist[z] = 1
    text = ''
    Z = hist.keys()
    Z.sort()
    for z in Z:
        text += symbols[z]
        n = hist[z]
        if n > 1:
            text += '<sub>%d</sub>' % n
    return text

class Status:
    def __init__(self, vbox):
        self.eventbox = gtk.EventBox()
        self.label = gtk.Label()
        self.eventbox.add(self.label)
        self.label.show()
        if gtk.pygtk_version < (2, 12):
            self.set_tip(self.eventbox, _('Tip for status box ...'))
        else:
            self.eventbox.set_tooltip_text(_('Tip for status box ...'))
        pack(vbox, self.eventbox)
        self.ordered_indices = []

    def status(self):
        # use where here:  XXX
        indices = np.arange(self.images.natoms)[self.images.selected]
        ordered_indices = self.images.selected_ordered
        n = len(indices)
        self.nselected = n
        
        if n == 0:
            self.label.set_text('')
            return

        Z = self.images.Z[indices]
        R = self.R[indices]

        if n == 1:
            tag = self.images.T[self.frame,indices][0]
            mom = self.images.M[self.frame][indices]
            text = (u' #%d %s (%s): %.3f Å, %.3f Å, %.3f Å ' %
                    ((indices[0], names[Z[0]], symbols[Z[0]]) + tuple(R[0])))
            # TRANSLATORS: mom refers to magnetic moment
            text += _(' tag=%(tag)s mom=%(mom)1.2f') % dict(tag=tag, mom=mom)
        elif n == 2:
            D = R[0] - R[1]
            d = sqrt(np.dot(D, D))
            text = u' %s-%s: %.3f Å' % (symbols[Z[0]], symbols[Z[1]], d)
        elif n == 3:
            d = []
            for c in range(3):
                D = R[c] - R[(c + 1) % 3]
                d.append(np.dot(D, D))
            a = []
            for c in range(3):
                t1 = 0.5 * (d[c] + d[(c + 1) % 3] - d[(c + 2) % 3])
                t2 = sqrt(d[c] * d[(c + 1) % 3])
                try:
                    t3 = acos(t1 / t2)
                except ValueError:
                    if t1 > 0:
                        t3 = 0
                    else:
                        t3 = pi
                a.append(t3 * 180 / pi)
            text = (u' %s-%s-%s: %.1f°, %.1f°, %.1f°' %
                    tuple([symbols[z] for z in Z] + a))
        elif len(ordered_indices) == 4:
            R = self.R[ordered_indices]
            Z = self.images.Z[ordered_indices]
            a    = R[1]-R[0]
            b    = R[2]-R[1]
            c    = R[3]-R[2]
            bxa  = np.cross(b,a)
            bxa /= np.sqrt(np.vdot(bxa,bxa))
            cxb  = np.cross(c,b)
            cxb /= np.sqrt(np.vdot(cxb,cxb))
            angle = np.vdot(bxa,cxb)
            if angle < -1: angle = -1
            if angle >  1: angle =  1
            angle = np.arccos(angle)
            if (np.vdot(bxa,c)) > 0: angle = 2*np.pi-angle
            angle = angle*180.0/np.pi
            text = (u'%s %s->%s->%s->%s: %.1f°'
                    % tuple([_('dihedral')] + [symbols[z] for z in Z]+[angle]))
        else:
            text = ' ' + formula(Z)
            
        self.label.set_markup(text)
        
if __name__ == '__main__':
    import os
    os.system('python gui.py')
