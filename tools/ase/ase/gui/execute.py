#!/usr/bin/env python
import __future__
import gtk
from gettext import gettext as _
import os.path
import numpy as np
import sys

from ase.gui.widgets import pack, Help
from ase.data.colors import jmol_colors
from ase.atoms import Atoms

class Execute(gtk.Window):
    """ The Execute class provides an expert-user window for modification
    and evaluation of system properties with a simple one-line command structure.
    There are two types of commands, one set only applies to the global image and
    one set applies to all atoms. If the command line contains any of the atom
    commands, then it is executed separately for all atoms and for all images.
    Otherwise it is executed only once per image. 

    Please do not mix global and atom commands."""
    
    terminal_help_txt=_("""
    Global commands work on all frames or only on the current frame
    - Assignment of a global variable may not reference a local one
    - use 'Current frame' switch to switch off application to all frames
    <c>e</c>:\t\ttotal energy of one frame
    <c>fmax</c>:\tmaximal force in one frame
    <c>A</c>:\tunit cell
    <c>E</c>:\t\ttotal energy array of all frames
    <c>F</c>:\t\tall forces in one frame
    <c>M</c>:\tall magnetic moments
    <c>R</c>:\t\tall atomic positions
    <c>S</c>:\tall selected atoms (boolean array)
    <c>D</c>:\tall dynamic atoms (boolean array)
    examples: <c>frame = 1</c>, <c>A[0][1] += 4</c>, <c>e-E[-1]</c>

    Atom commands work on each atom (or a selection) individually
    - these can use global commands on the RHS of an equation
    - use 'selected atoms only' to restrict application of command
    <c>x,y,z</c>:\tatomic coordinates
    <c>r,g,b</c>:\tatom display color, range is [0..1]
    <c>rad</c>:\tatomic radius for display
    <c>s</c>:\t\tatom is selected
    <c>d</c>:\t\tatom is movable
    <c>f</c>:\t\tforce
    <c>Z</c>:\tatomic number
    <c>m</c>:\tmagnetic moment
    examples: <c>x -= A[0][0], s = z > 5, Z = 6</c>

    Special commands and objects:
    <c>sa,cf</c>:\t(un)restrict to selected atoms/current frame
    <c>frame</c>:\tframe number
    <c>center</c>:\tcenters the system in its existing unit cell
    <c>del S</c>:\tdelete selection
    <c>CM</c>:\tcenter of mass
    <c>ans[-i]</c>:\tith last calculated result
    <c>exec file</c>: executes commands listed in file
    <c>cov[Z]</c>:(read only): covalent radius of atomic number Z
    <c>gui</c>:\tadvanced: ag window python object
    <c>img</c>:\tadvanced: ag images object
    """)
    
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.gui = gui
        self.set_title(_('Expert user mode'))
        vbox = gtk.VBox()
        vbox.set_border_width(5)
        self.sw = gtk.ScrolledWindow()
        self.sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.textview = gtk.TextView()
        self.textbuffer = self.textview.get_buffer()
        self.textview.set_editable(False)
        self.textview.set_cursor_visible(False)
        self.sw.add(self.textview)
        pack(vbox, self.sw, expand=True, padding = 5)
        self.sw.set_size_request(540, 150)
        self.textview.show()
        self.add_text(_('Welcome to the ASE Expert user mode'))
        self.cmd = gtk.Entry(60)
        self.cmd.connect('activate', self.execute)
        self.cmd.connect('key-press-event', self.update_command_buffer)
        pack(vbox, [gtk.Label('>>>'),self.cmd])
        self.cmd_buffer = getattr(gui,'expert_mode_buffer',[''])
        self.cmd_position = len(self.cmd_buffer)-1
        self.selected = gtk.CheckButton(_('Only selected atoms (sa)   '))
        self.selected.connect('toggled',self.selected_changed)
        self.images_only = gtk.CheckButton(_('Only current frame (cf)  '))
        self.images_only.connect('toggled',self.images_changed)
        pack(vbox, [self.selected, self.images_only])
        save_button = gtk.Button(stock=gtk.STOCK_SAVE)
        save_button.connect('clicked',self.save_output)
        help_button = gtk.Button(stock=gtk.STOCK_HELP)
        help_button.connect('clicked',self.terminal_help,"")
        stop_button = gtk.Button(stock=gtk.STOCK_STOP)
        stop_button.connect('clicked',self.stop_execution)
        self.stop = False
        pack(vbox, [gtk.Label(_('Global: Use A, D, E, M, N, R, S, n, frame;'
                                ' Atoms: Use a, f, m, s, x, y, z, Z     ')),
                    stop_button, help_button, save_button], end = True)
        self.add(vbox)
        vbox.show()
        self.show()
        # set color mode to manual when opening this window for rgb manipulation
        self.colors = self.gui.get_colors()
        rgb_data = self.gui.get_colors(rgb = True)
        self.rgb_data = []  # ensure proper format of rgb_data
        for i, rgb in enumerate(rgb_data):
            self.rgb_data += [[i, rgb]]
        self.gui.colordata = self.rgb_data
        self.gui.colors = list(self.colors)
        self.gui.colormode = 'manual'        
        self.cmd.grab_focus()

    def execute(self, widget=None, cmd = None):
        global_commands = ['A','Col','D','e','E','F','frame','M','n','N','R','S']  # explicitly 'implemented' commands for use on whole system or entire single frame
        index_commands  = ['a','b','d','f','g','m','r','rad','s','x','y','z','Z']  # commands for use on all (possibly selected) atoms

        new = self.gui.drawing_area.window.new_gc
        alloc = self.gui.colormap.alloc_color

        self.stop = False
        if cmd is None:
            cmd = self.cmd.get_text().strip()
            if len(cmd) == 0:
                return
            self.add_text('>>> '+cmd)
            self.cmd_buffer[-1] = cmd
            self.cmd_buffer += ['']
            setattr(self.gui,'expert_mode_buffer', self.cmd_buffer)
            self.cmd_position = len(self.cmd_buffer)-1
            self.cmd.set_text('')
        else:
            self.add_text('--> '+cmd)

        gui = self.gui
        img = gui.images
        frame = gui.frame
        N = img.nimages
        n = img.natoms
        S = img.selected
        D = img.dynamic[:, np.newaxis]
        E = img.E
        if self.selected.get_active():
            indices = np.where(S)[0]
        else:
            indices = range(n)

        ans = getattr(gui,'expert_mode_answers',[])

        loop_images = range(N)
        if self.images_only.get_active():
            loop_images = [self.gui.frame]

        # split off the first valid command in cmd to determine whether
        # it is global or index based, this includes things such as 4*z and z*4
        index_based = False
        first_command = cmd.split()[0]
        special = ['=',',','+','-','/','*',';','.','[',']','(',')',
                   '{','}','0','1','2','3','4','5','6','7','8','9']
        while first_command[0] in special and len(first_command)>1:
            first_command = first_command[1:]
        for c in special:
            if c in first_command:
                first_command = first_command[:first_command.find(c)]
        for c in index_commands:
            if c == first_command:
                index_based = True

        name = os.path.expanduser('~/.ase/'+cmd)
        # check various special commands: 
        if os.path.exists(name):   # run script from default directory
            self.run_script(name)
        elif cmd == 'del S':       # delete selection
            gui.delete_selected_atoms()
        elif cmd == 'sa':          # selected atoms only
            self.selected.set_active(not self.selected.get_active())
        elif cmd == 'cf':          # current frame only
            self.images_only.set_active(not self.images_only.get_active())
        elif cmd == 'center':      # center system
            img.center()
        elif cmd == 'CM':          # calculate center of mass
            for i in loop_images:
                if self.stop:
                    break
                atoms = Atoms(positions=img.P[i][indices],
                              numbers=img.Z[indices])
                self.add_text(repr(atoms.get_center_of_mass()))
                ans += [atoms.get_center_of_mass()]
        elif first_command == 'exec': # execute script
            name = cmd.split()[1]
            if '~' in name:
                name = os.path.expanduser(name)
            if os.path.exists(name):
                self.run_script(name)
            else:
                self.add_text(_('*** WARNING: file does not exist - %s') % name)
        else:
            code = compile(cmd + '\n', 'execute.py', 'single',
                           __future__.CO_FUTURE_DIVISION)
            if index_based and len(indices) == 0 and self.selected.get_active():
                self.add_text(_("*** WARNING: No atoms selected to work with"))
            for i in loop_images:
                if self.stop:
                    break
                R = img.P[i][indices]
                A = img.A[i]
                F = img.F[i][indices]
                e = img.E[i]
                M = img.M[i][indices]
                Col = []
                cov = img.covalent_radii
                for j in indices:
                    Col += [gui.colordata[j]]
                if len(indices) > 0:
                    fmax = max(((F * D[indices])**2).sum(1)**.5)
                else:
                    fmax = None
                frame = gui.frame
            
                if not index_based:
                    try:
                        self.add_text(repr(eval(cmd)))
                        ans += [eval(cmd)]
                    except:
                        exec code
                    gui.set_frame(frame)
                    if gui.movie_window is not None:
                        gui.movie_window.frame_number.value = frame
                    img.selected      = S
                    img.A[i]          = A
                    img.P[i][indices] = R
                    img.M[i][indices] = M
                else:
                    for n,a in enumerate(indices):
                        if self.stop:
                            break
                        x, y, z = R[n]
                        r, g, b = Col[n][1]
                        d = D[a]
                        f = np.vdot(F[n]*d,F[n]*d)**0.5
                        s = S[a]
                        Z = img.Z[a]
                        Zold = Z
                        m = M[n]
                        rad = img.r[a]
                        try:
                            self.add_text(repr(eval(cmd)))
                            ans += [eval(cmd)]
                        except:
                            exec code
                        S[a] = s
                        img.P[i][a] = x, y, z
                        img.Z[a] = Z
                        img.r[a] = rad
                        img.dynamic[a] = d
                        if Z != Zold:
                            img.r[a] = cov[Z] * 0.89
                            r,g,b = jmol_colors[Z]
                        gui.colordata[a] = [a,[r,g,b]]                            
                        color = tuple([int(65535*x) for x in [r,g,b]])
                        gui.colors[a] = new(alloc(*color))
                        img.M[i][a] = m
        setattr(self.gui,'expert_mode_answers', ans)
        gui.set_frame(frame,init=True)

    def add_text(self,val):
        text_end = self.textbuffer.get_end_iter()
        self.textbuffer.insert(text_end,val+'\n');
        if self.sw.get_vscrollbar() is not None:
            scroll = self.sw.get_vscrollbar().get_adjustment()
            scroll.set_value(scroll.get_upper())
        
    def selected_changed(self, *args):
        if self.selected.get_active():
            self.add_text(_('*** Only working on selected atoms'))
        else:
            self.add_text(_('*** Working on all atoms'))

    def images_changed(self, *args):
        if self.images_only.get_active():
            self.add_text(_('*** Only working on current image'))
        else:
            self.add_text(_('*** Working on all images'))

    def update_command_buffer(self, entry, event, *args):
        arrow = {gtk.keysyms.Up: -1, gtk.keysyms.Down: 1}.get(event.keyval, None)
        if arrow is not None:
            self.cmd_position += arrow
            self.cmd_position = max(self.cmd_position,0)
            self.cmd_position = min(self.cmd_position,len(self.cmd_buffer)-1)
            cmd = self.cmd_buffer[self.cmd_position]
            self.cmd.set_text(cmd)
            return True
        else:
            return False

    def save_output(self, *args):
        chooser = gtk.FileChooserDialog(
            _('Save Terminal text ...'), None, gtk.FILE_CHOOSER_ACTION_SAVE,
            (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
             gtk.STOCK_SAVE, gtk.RESPONSE_OK))
        save = chooser.run()
        if save == gtk.RESPONSE_OK or save == gtk.RESPONSE_SAVE:
            filename = chooser.get_filename()
            text = self.textbuffer.get_text(self.textbuffer.get_start_iter(),
                                            self.textbuffer.get_end_iter())
            fd = open(filename,'w')
            fd.write(text)
            fd.close()
            chooser.destroy()

    def run_script(self, name):
        commands = open(name,'r').readlines()
        for c_parse in commands:
            c = c_parse.strip()
            if '#' in c:
                c = c[:c.find('#')].strip()
            if len(c) > 0:
                self.execute(cmd = c.strip())
            
    def terminal_help(self,*args):
        Help(self.terminal_help_txt)

    def stop_execution(self, *args):
        self.stop = True
        
    python = execute
