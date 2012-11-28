# encoding: utf-8
"""colors.py - select how to color the atoms in the GUI."""


import gtk
from gettext import gettext as _
from ase.gui.widgets import pack, cancel_apply_ok, oops, help
import ase
from ase.data.colors import jmol_colors
import numpy as np
import colorsys

named_colors = ('Green', 'Yellow', 'Blue', 'Red', 'Orange', 'Cyan',
                'Magenta', 'Black', 'White', 'Grey', 'Violet', 'Brown',
                'Navy')

class ColorWindow(gtk.Window):
    "A window for selecting how to color the atoms."
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.gui = gui
        self.colormode = gui.colormode
        self.actual_colordata = None
        self.set_title(_("Colors"))
        vbox = gtk.VBox()
        self.add(vbox)
        vbox.show()
        # The main layout consists of two columns, the leftmost split in an upper and lower part.
        self.maintable = gtk.Table(2,2)
        pack(vbox, self.maintable)
        self.methodbox = gtk.VBox()
        self.methodbox.show()
        self.maintable.attach(self.methodbox, 0, 1, 0, 1)
        self.scalebox = gtk.VBox()
        self.scalebox.show()
        self.maintable.attach(self.scalebox, 0, 1, 1, 2)
        self.colorbox = gtk.Frame()
        self.colorbox.show()
        self.maintable.attach(self.colorbox, 1, 2, 0, 2, gtk.EXPAND)
        # Upper left: Choose how the atoms are colored.
        lbl = gtk.Label(_("Choose how the atoms are colored:"))
        pack(self.methodbox, [lbl])
        self.radio_jmol = gtk.RadioButton(None, _('By atomic number, default "jmol" colors'))
        self.radio_atno = gtk.RadioButton(self.radio_jmol,
                                          _('By atomic number, user specified'))
        self.radio_tag = gtk.RadioButton(self.radio_jmol, _('By tag'))
        self.radio_force = gtk.RadioButton(self.radio_jmol, _('By force'))
        self.radio_velocity = gtk.RadioButton(self.radio_jmol, _('By velocity'))
        self.radio_manual = gtk.RadioButton(self.radio_jmol, _('Manually specified'))
        self.radio_same = gtk.RadioButton(self.radio_jmol, _('All the same color'))
        self.force_box = gtk.VBox()
        self.velocity_box = gtk.VBox()
        for widget in (self.radio_jmol, self.radio_atno, self.radio_tag,
                      self.radio_force, self.force_box, self.radio_velocity,
                      self.velocity_box, self.radio_manual, self.radio_same):
            pack(self.methodbox, [widget])
            if isinstance(widget, gtk.RadioButton):
                widget.connect('toggled', self.method_radio_changed)
        # Now fill in the box for additional information in case the force is used.
        self.force_label = gtk.Label(_("This should not be displayed!"))
        pack(self.force_box, [self.force_label])
        self.force_min = gtk.Adjustment(0.0, 0.0, 100.0, 0.05)
        self.force_max = gtk.Adjustment(0.0, 0.0, 100.0, 0.05)
        self.force_steps = gtk.Adjustment(10, 2, 500, 1)
        force_apply = gtk.Button(_('Update'))
        force_apply.connect('clicked', self.set_force_colors)
        pack(self.force_box, [gtk.Label(_('Min: ')),
                              gtk.SpinButton(self.force_min, 10.0, 2),
                              gtk.Label(_('  Max: ')),
                              gtk.SpinButton(self.force_max, 10.0, 2),
                              gtk.Label(_('  Steps: ')),
                              gtk.SpinButton(self.force_steps, 1, 0),
                              gtk.Label('  '),
                              force_apply])
        self.force_box.hide()
        # Now fill in the box for additional information in case the velocity is used.
        self.velocity_label = gtk.Label("This should not be displayed!")
        pack(self.velocity_box, [self.velocity_label])
        self.velocity_min = gtk.Adjustment(0.0, 0.0, 10.0, 0.005)
        self.velocity_max = gtk.Adjustment(0.0, 0.0, 10.0, 0.005)
        self.velocity_steps = gtk.Adjustment(10, 2, 500, 1)
        velocity_apply = gtk.Button(_('Update'))
        velocity_apply.connect('clicked', self.set_velocity_colors)
        pack(self.velocity_box, [gtk.Label(_('Min: ')),
                                 gtk.SpinButton(self.velocity_min, 10.0, 3),
                                 gtk.Label(_('  Max: ')),
                                 gtk.SpinButton(self.velocity_max, 10.0, 3),
                                 gtk.Label(_('  Steps: ')),
                                 gtk.SpinButton(self.velocity_steps, 1, 0),
                                 gtk.Label('  '),
                                 velocity_apply])
        self.velocity_box.hide()
        # Lower left: Create a color scale
        pack(self.scalebox, gtk.Label(""))
        lbl = gtk.Label(_('Create a color scale:'))
        pack(self.scalebox, [lbl])
        color_scales = (
            _('Black - white'),
            _('Black - red - yellow - white'),
            _('Black - green - white'),
            _('Black - blue - cyan'),
            _('Hue'),
            _('Named colors')
            )
        self.scaletype_created = None
        self.scaletype = gtk.combo_box_new_text()
        for s in color_scales:
            self.scaletype.append_text(s)
        self.createscale = gtk.Button(_("Create"))
        pack(self.scalebox, [self.scaletype, self.createscale])
        self.createscale.connect('clicked', self.create_color_scale)
        # The actually colors are specified in a box possibly with scrollbars
        self.colorwin = gtk.ScrolledWindow()
        self.colorwin.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
        self.colorwin.show()
        self.colorbox.add(self.colorwin)
        self.colorwin.add_with_viewport(gtk.VBox()) # Dummy contents
        buts = cancel_apply_ok(cancel=lambda widget: self.destroy(),
                               apply=self.apply,
                               ok=self.ok)
        pack(vbox, [buts], end=True, bottom=True)
        # Make the initial setup of the colors
        self.color_errors = {}
        self.init_colors_from_gui()
        self.show()
        gui.register_vulnerable(self)

    def notify_atoms_changed(self):
        "Called by gui object when the atoms have changed."
        self.destroy()
  
    def init_colors_from_gui(self):
        cm = self.gui.colormode
        # Disallow methods if corresponding data is not available
        if not self.gui.images.T.any():
            self.radio_tag.set_sensitive(False)
            if self.radio_tag.get_active() or cm == 'tag':
                self.radio_jmol.set_active(True)
                return
        else:
            self.radio_tag.set_sensitive(True)
        if np.isnan(self.gui.images.F).any() or not self.gui.images.F.any():
            self.radio_force.set_sensitive(False)
            if self.radio_force.get_active() or cm == 'force':
                self.radio_jmol.set_active(True)
                return
        else:
            self.radio_force.set_sensitive(True)
        if np.isnan(self.gui.images.V).any() or not self.gui.images.V.any():
            self.radio_velocity.set_sensitive(False)
            if self.radio_velocity.get_active() or cm == 'velocity':
                self.radio_jmol.set_active(True)
                return
        else:
            self.radio_velocity.set_sensitive(True)
        self.radio_manual.set_sensitive(self.gui.images.natoms <= 1000)
        # Now check what the current color mode is
        if cm == 'jmol':
            self.radio_jmol.set_active(True)
            self.set_jmol_colors()
        elif cm == 'Z':
            self.radio_atno.set_active(True)
        elif cm == 'tag':
            self.radio_tag.set_active(True)
        elif cm == 'force':
            self.radio_force.set_active(True)
        elif cm == 'velocity':
            self.radio_velocity.set_active(True)
        elif cm == 'manual':
            self.radio_manual.set_active(True)
        elif cm == 'same':
            self.radio_same.set_active(True)
            
    def method_radio_changed(self, widget=None):
        "Called when a radio button is changed."
        self.scaletype_created = None
        self.scaletype.set_active(-1)
        if not widget.get_active():
            # Ignore most events when a button is turned off.
            if widget is self.radio_force:
                self.force_box.hide()
            if widget is self.radio_velocity:
                self.velocity_box.hide()
            return  
        if widget is self.radio_jmol:
            self.set_jmol_colors()
        elif widget is self.radio_atno:
            self.set_atno_colors()
        elif widget is self.radio_tag:
            self.set_tag_colors()
        elif widget is self.radio_force:
            self.show_force_stuff()
            self.set_force_colors()
        elif widget is self.radio_velocity:
            self.show_velocity_stuff()
            self.set_velocity_colors()
        elif widget is self.radio_manual:
            self.set_manual_colors()
        elif widget is self.radio_same:
            self.set_same_color()
        else:
            raise RuntimeError('Unknown widget in method_radio_changed')
            
    def make_jmol_colors(self):
        "Set the colors to the default jmol colors"
        self.colordata_z = []
        hasfound = {}
        for z in self.gui.images.Z:
            if z not in hasfound:
                hasfound[z] = True
                self.colordata_z.append([z, jmol_colors[z]])

    def set_jmol_colors(self):
        "We use the immutable jmol colors."
        self.make_jmol_colors()
        self.set_atno_colors()
        for entry in self.color_entries:
            entry.set_sensitive(False)
        self.colormode = 'jmol'
        
    def set_atno_colors(self):
        "We use user-specified per-element colors."
        if not hasattr(self, 'colordata_z'):
            # No initial colors.  Use jmol colors
            self.make_jmol_colors()
        self.actual_colordata = self.colordata_z
        self.color_labels = ["%i (%s):" % (z, ase.data.chemical_symbols[z])
                             for z, col in self.colordata_z]
        self.make_colorwin()
        self.colormode = "atno"

    def set_tag_colors(self):
        "We use per-tag colors."
        # Find which tags are in use
        tags = self.gui.images.T
        existingtags = []
        for t in range(tags.min(), tags.max()+1):
            if t in tags:
                existingtags.append(t)
        if not hasattr(self, 'colordata_tags') or len(self.colordata_tags) != len(existingtags):
            colors = self.get_named_colors(len(existingtags))
            self.colordata_tags = [[x, y] for x, y in
                                   zip(existingtags, colors)]
        self.actual_colordata = self.colordata_tags
        self.color_labels = [str(x)+':' for x, y in self.colordata_tags]
        self.make_colorwin()
        self.colormode = 'tags'

    def set_same_color(self):
        "All atoms have the same color"
        if not hasattr(self, 'colordata_same'):
            try:
                self.colordata_same = self.actual_colordata[0:1]
            except AttributeError:
                self.colordata_same = self.get_named_colors(1)
        self.actual_colordata = self.colordata_same
        self.actual_colordata[0][0] = 0
        self.color_labels = ['all:']
        self.make_colorwin()
        self.colormode = 'same'

    def set_force_colors(self, *args):
        "Use the forces as basis for the colors."
        borders = np.linspace(self.force_min.value,
                              self.force_max.value,
                              self.force_steps.value,
                              endpoint=False)
        if self.scaletype_created is None:
            colors = self.new_color_scale([[0, [1,1,1]],
                                           [1, [0,0,1]]], len(borders))
        elif (not hasattr(self, 'colordata_force') or
            len(self.colordata_force) != len(borders)):
            colors = self.get_color_scale(len(borders), self.scaletype_created)
        else:
            colors = [y for x, y in self.colordata_force]
        self.colordata_force = [[x, y] for x, y in zip(borders, colors)]
        self.actual_colordata = self.colordata_force
        self.color_labels = ["%.2f:" % x for x, y in self.colordata_force]
        self.make_colorwin()
        self.colormode = 'force'
        fmin = self.force_min.value
        fmax = self.force_max.value
        factor = self.force_steps.value / (fmax -fmin)
        self.colormode_force_data = (fmin, factor)

    def set_velocity_colors(self, *args):
        "Use the velocities as basis for the colors."
        borders = np.linspace(self.velocity_min.value,
                              self.velocity_max.value,
                              self.velocity_steps.value,
                              endpoint=False)
        if self.scaletype_created is None:
            colors = self.new_color_scale([[0, [1,1,1]],
                                           [1, [1,0,0]]], len(borders))
        elif (not hasattr(self, 'colordata_velocity') or
            len(self.colordata_velocity) != len(borders)):
            colors = self.get_color_scale(len(borders), self.scaletype_created)
        else:
            colors = [y for x, y in self.colordata_velocity]
        self.colordata_velocity = [[x, y] for x, y in zip(borders, colors)]
        self.actual_colordata = self.colordata_velocity
        self.color_labels = ["%.2f:" % x for x, y in self.colordata_velocity]
        self.make_colorwin()
        self.colormode = 'velocity'
        vmin = self.velocity_min.value
        vmax = self.velocity_max.value
        factor = self.velocity_steps.value / (vmax -vmin)
        self.colormode_velocity_data = (vmin, factor)

    def set_manual_colors(self):
        "Set colors of all atoms from the last selection."
        # We cannot directly make np.arrays of the colors, as they may
        # be sequences of the same length, causing creation of a 2D
        # array of characters/numbers instead of a 1D array of
        # objects.
        colors = np.array([None] * self.gui.images.natoms)
        if self.colormode in ['atno', 'jmol', 'tags']:
            maxval = max([x for x, y in self.actual_colordata])
            oldcolors = np.array([None] * (maxval+1))
            for x, y in self.actual_colordata:
                oldcolors[x] = y
            if self.colormode == 'tags':
                colors[:] = oldcolors[self.gui.images.T[self.gui.frame]]
            else:
                colors[:] = oldcolors[self.gui.images.Z]
        elif self.colormode == 'force':
            oldcolors = np.array([None] * len(self.actual_colordata))
            oldcolors[:] = [y for x, y in self.actual_colordata]
            F = self.gui.images.F[self.gui.frame]
            F = np.sqrt((F * F).sum(axis=-1))
            nF = (F - self.colormode_force_data[0]) * self.colormode_force_data[1]
            nF = np.clip(nF.astype(int), 0, len(oldcolors)-1)
            colors[:] = oldcolors[nF]
        elif self.colormode == 'velocity':
            oldcolors = np.array([None] * len(self.actual_colordata))
            oldcolors[:] = [y for x, y in self.actual_colordata]
            V = self.gui.images.V[self.gui.frame]
            V = np.sqrt((V * V).sum(axis=-1))
            nV = (V - self.colormode_velocity_data[0]) * self.colormode_velocity_data[1]
            nV = np.clip(nV.astype(int), 0, len(oldcolors)-1)
            colors[:] = oldcolors[nV]
        elif self.colormode == 'same':
            oldcolor = self.actual_colordata[0][1]
            if len(colors) == len(oldcolor):
                # Direct assignment would be e.g. one letter per atom. :-(
                colors[:] = [oldcolor] * len(colors)
            else:
                colors[:] = oldcolor
        elif self.colormode == 'manual':
            if self.actual_colordata is None:   # import colors from gui, if they don't exist already
                colors = [y for x,y in self.gui.colordata]

        self.color_labels = ["%d:" % i for i in range(len(colors))]
        self.actual_colordata = [[i, x] for i, x in enumerate(colors)]
        self.make_colorwin()
        self.colormode = 'manual'

    def show_force_stuff(self):
        "Show and update widgets needed for selecting the force scale."
        self.force_box.show()
        F = np.sqrt(((self.gui.images.F*self.gui.images.dynamic[:,np.newaxis])**2).sum(axis=-1))
        fmax = F.max()
        nimages = self.gui.images.nimages
        assert len(F) == nimages
        if nimages > 1:
            fmax_frame = self.gui.images.F[self.gui.frame].max()
            txt = _("Max force: %.2f (this frame), %.2f (all frames)") % (fmax_frame, fmax)
        else:
            txt = _("Max force: %.2f.") % (fmax,)
        self.force_label.set_text(txt)
        if self.force_max.value == 0.0:
            self.force_max.value = fmax

    def show_velocity_stuff(self):
        "Show and update widgets needed for selecting the velocity scale."
        self.velocity_box.show()
        V = np.sqrt((self.gui.images.V * self.gui.images.V).sum(axis=-1))
        vmax = V.max()
        nimages = self.gui.images.nimages
        assert len(V) == nimages
        if nimages > 1:
            vmax_frame = self.gui.images.V[self.gui.frame].max()
            txt = _("Max velocity: %.2f (this frame), %.2f (all frames)") % (vmax_frame, vmax)
        else:
            txt = _("Max velocity: %.2f.") % (vmax,)
        self.velocity_label.set_text(txt)
        if self.velocity_max.value == 0.0:
            self.velocity_max.value = vmax
        
    def make_colorwin(self):
        """Make the list of editable color entries.

        Uses self.actual_colordata and self.color_labels.  Produces self.color_entries.
        """
        assert len(self.actual_colordata) == len(self.color_labels)
        self.color_entries = []
        old = self.colorwin.get_child()
        self.colorwin.remove(old)
        del old
        table = gtk.Table(len(self.actual_colordata)+1, 4)
        self.colorwin.add_with_viewport(table)
        table.show()
        self.color_display = []
        for i in range(len(self.actual_colordata)):
            lbl = gtk.Label(self.color_labels[i])
            entry = gtk.Entry(max=20)
            val = self.actual_colordata[i][1]
            error = False
            if not isinstance(val, str):
                assert len(val) == 3
                intval = tuple(np.round(65535*np.array(val)).astype(int))
                val = "%.3f, %.3f, %.3f" % tuple(val)
                clr = gtk.gdk.Color(*intval)
            else:
                try:
                    clr = gtk.gdk.color_parse(val)
                except ValueError:
                    error = True
            entry.set_text(val)
            blob = gtk.EventBox()
            space = gtk.Label
            space = gtk.Label("    ")
            space.show()
            blob.add(space)
            if error:
                space.set_text(_("ERROR"))
            else:
                blob.modify_bg(gtk.STATE_NORMAL, clr)
            table.attach(lbl, 0, 1, i, i+1, yoptions=0)
            table.attach(entry, 1, 2, i, i+1, yoptions=0)
            table.attach(blob, 2, 3, i, i+1, yoptions=0)
            lbl.show()
            entry.show()
            blob.show()
            entry.connect('activate', self.entry_changed, i)
            self.color_display.append(blob)
            self.color_entries.append(entry)
            
    def entry_changed(self, widget, index):
        """The user has changed a color."""
        txt = widget.get_text()
        txtfields = txt.split(',')
        if len(txtfields) == 3:
            self.actual_colordata[index][1] = [float(x) for x in txtfields]
            val = tuple([int(65535*float(x)) for x in txtfields])
            clr = gtk.gdk.Color(*val)
        else:
            self.actual_colordata[index][1] = txt
            try:
                clr = gtk.gdk.color_parse(txt)
            except ValueError:
                # Cannot parse the color
                displ = self.color_display[index]
                displ.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse('white'))
                displ.get_child().set_text(_("ERR"))
                self.color_errors[index] = (self.color_labels[index], txt)
                return
        self.color_display[index].get_child().set_text("    ") # Clear error message
        self.color_errors.pop(index, None)
        self.color_display[index].modify_bg(gtk.STATE_NORMAL, clr)
        
    def create_color_scale(self, *args):
        if self.radio_jmol.get_active():
            self.radio_atno.set_active(1)
        n = len(self.color_entries)
        s = self.scaletype.get_active()
        scale = self.get_color_scale(n, s)
        self.scaletype_created = s
        for i in range(n):
            if isinstance(scale[i], str):
                self.color_entries[i].set_text(scale[i])
            else:
                s = "%.3f, %.3f, %.3f" % tuple(scale[i])
                self.color_entries[i].set_text(s)
            self.color_entries[i].activate()

    def get_color_scale(self, n, s):
        if s == 0:
            # Black - White
            scale = self.new_color_scale([[0, [0,0,0]],
                                          [1, [1,1,1]]], n)
        elif s == 1:
            # Black - Red - Yellow - White (STM colors)
            scale = self.new_color_scale([[0, [0,0,0]],
                                          [0.33, [1,0,0]],
                                          [0.67, [1,1,0]],
                                          [1, [1,1,1]]], n)
        elif s == 2:
            # Black - Green - White
            scale = self.new_color_scale([[0, [0,0,0]],
                                          [0.5, [0,0.9,0]],
                                          [0.75, [0.2,1.0,0.2]],
                                          [1, [1,1,1]]], n)
        elif s == 3:
            # Black - Blue - Cyan
            scale = self.new_color_scale([[0, [0,0,0]],
                                          [0.5, [0,0,1]],
                                          [1, [0,1,1]]], n)
        elif s == 4:
            # Hues
            hues = np.linspace(0.0, 1.0, n, endpoint=False)
            scale = ["%.3f, %.3f, %.3f" % colorsys.hls_to_rgb(h, 0.5, 1)
                     for h in hues]
        elif s == 5:
            # Named colors
            scale = self.get_named_colors(n)
        else:
            scale = None
        return scale

    def new_color_scale(self, fixpoints, n):
        "Create a homogeneous color scale."
        x = np.array([a[0] for a in fixpoints], float)
        y = np.array([a[1] for a in fixpoints], float)
        assert y.shape[1] == 3
        res = []
        for a in np.linspace(0.0, 1.0, n, endpoint=True):
            n = x.searchsorted(a)
            if n == 0:
                v = y[0]  # Before the start
            elif n == len(x):
                v = x[-1] # After the end
            else:
                x0 = x[n-1]
                x1 = x[n]
                y0 = y[n-1]
                y1 = y[n]
                v = y0 + (y1 - y0) / (x1 - x0) * (a - x0)
            res.append(v)
        return res

    def get_named_colors(self, n):
        if n <= len(named_colors):
            return named_colors[:n]
        else:
            return named_colors + ('Black',) * (n - len(named_colors))
        
    def apply(self, *args):
        #if self.colormode in ['atno', 'jmol', 'tags']:
        # Color atoms according to an integer value number
        if self.color_errors:
            oops(_("Incorrect color specification"),
                 "%s: %s" % self.color_errors.values()[0])
            return False
        colordata = self.actual_colordata
        if self.colormode == 'force':
            # Use integers instead for border values
            colordata = [[i, x[1]] for i, x in enumerate(self.actual_colordata)]
            self.gui.colormode_force_data = self.colormode_force_data
        if self.colormode == 'velocity':
            # Use integers instead for border values
            colordata = [[i, x[1]] for i, x in enumerate(self.actual_colordata)]
            self.gui.colormode_velocity_data = self.colormode_velocity_data
        maxval = max([x for x, y in colordata])
        self.gui.colors = [None] * (maxval + 1)
        new = self.gui.drawing_area.window.new_gc
        alloc = self.gui.colormap.alloc_color
        for z, val in colordata:
            if isinstance(val, str):
                self.gui.colors[z] = new(alloc(val))
            else:
                clr = tuple([int(65535*x) for x in val])
                assert len(clr) == 3
                self.gui.colors[z] = new(alloc(*clr))
        self.gui.colormode = self.colormode
        self.gui.colordata = self.actual_colordata
        self.gui.draw()
        return True

    def cancel(self, *args):
        self.destroy()

    def ok(self, *args):
        if self.apply():
            self.destroy()
        
