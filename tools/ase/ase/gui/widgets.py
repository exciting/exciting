from gettext import gettext as _
import gtk
import re

class Number(gtk.SpinButton):
    def __init__(self, value=0,
                 lower=0, upper=10000,
                 step_incr=1, page_incr=10,
                 climb_rate=0.5, digits=0):
        self.adj = gtk.Adjustment(value, lower, upper, step_incr, page_incr, 0)
        gtk.SpinButton.__init__(self, self.adj, climb_rate, digits)

    def connect(self, *args):
        return self.adj.connect(*args)


class Menu:
    def __init__(self, menubar, name, items):
        self.items = {}
        menu = gtk.Menu()
        for data in items:
            text = data[0]
            callback = data[1]
            args = data[2:]
            menuitem = gtk.MenuItem(text)
            menu.append(menuitem)
            menuitem.connect('activate', callback, *args)
            menuitem.show()
            self.items[text] = menuitem
        menuitem = gtk.MenuItem(name)
        menubar.append(menuitem)
        menuitem.set_submenu(menu)
        menuitem.show()


class Help(gtk.Window):
    __instance = None
    def __new__(cls, *args, **kwargs):
        # Make this a singleton.
        if Help.__instance is None:
            Help.__instance = gtk.Window.__new__(cls, *args, **kwargs)
        return Help.__instance
    def __init__(self, text):
        # Now, __init__ may be called multiple times!
        if not hasattr(self, '_initialized'):
            self.initialize(text)
        else:
            self.set_text(text)
        self.present()  # Show the window.
        
    def initialize(self, text):
        gtk.Window.__init__(self)
        self.set_title(_("Help"))
        self._initialized = True
        vbox = gtk.VBox()
        self.add(vbox)
        self.label = pack(vbox, gtk.Label())
        self.label.set_line_wrap(True)
        self.set_text(text)
        close = gtk.Button(_('Close'))
        pack(vbox, [close])
        close.connect('clicked', self.destroy)
        self.connect("delete-event", self.destroy) 
        self.show_all()

    def set_text(self, text):
        # Count line length without all the markup tags
        plaintext = ''.join(re.split('<[^>]+>', text))
        linelen = max([len(x) for x in plaintext.split('\n')])
        text = text.replace('<c>', '<span foreground="blue">')
        text = text.replace('</c>', '</span>')
        self.label.set_width_chars(linelen)
        self.label.set_line_wrap(False)
        self.label.set_markup(text)

    def destroy(self, *args):
        self.hide()
        return True  # Prevents destruction of the window.
        
def help(text):
    button = gtk.Button(_('Help'))
    button.connect('clicked', lambda widget, text=text: Help(text))
    return button


class Window(gtk.Window):
    def __init__(self, gui):
        self.gui = gui
        gtk.Window.__init__(self)
        self.set_title(_('Constraints'))
        vbox = gtk.VBox()
        b = pack(vbox, [gtk.Button(_('Constrain')),
                        gtk.Label(_(' selected atoms'))])[0]
        b.connect('clicked', self.selected)
        b = pack(vbox, [gtk.Button(_('Constrain')),
                        gtk.Label(_(' immobile atoms:'))])[0]
        b.connect('clicked', self.immobile)
        b = pack(vbox, gtk.Button(_('Clear constraint')))
        b.connect('clicked', self.clear)
        close = pack(vbox, gtk.Button(_('Close')))
        close.connect('clicked', lambda widget: self.destroy())
        self.add(vbox)
        vbox.show()
        self.show()

def pack(vbox, widgets, end=False, bottom=False, expand=False, padding=0):
    if not isinstance(widgets, list):
        widgets.show()
        if bottom:
            vbox.pack_end(widgets, expand, expand, padding)
        else:
            vbox.pack_start(widgets, expand, expand, padding)
        return widgets
    hbox = gtk.HBox(0, 0)
    hbox.show()
    if bottom:
        vbox.pack_end(hbox, expand, expand, padding)
    else:
        vbox.pack_start(hbox, expand, expand, padding)
    for widget in widgets:
        if type(widget) is gtk.Entry:
            widget.set_size_request(widget.get_max_length() * 9, 24)
        widget.show()
        if end and widget is widgets[-1]:
            hbox.pack_end(widget, expand, expand, padding)
        else:
            hbox.pack_start(widget, expand, expand, padding)
    return widgets

class cancel_apply_ok(gtk.HButtonBox):
    "Widget with Cancel, Apply and OK buttons.  The arguments are callbacks."
    def __init__(self, cancel, apply, ok):
        gtk.HButtonBox.__init__(self)
        cancel_but = gtk.Button(stock=gtk.STOCK_CANCEL)
        cancel_but.connect('clicked', cancel)
        apply_but = gtk.Button(stock=gtk.STOCK_APPLY)
        apply_but.connect('clicked', apply)
        ok_but = gtk.Button(stock=gtk.STOCK_OK)
        ok_but.connect('clicked', ok)
        for w in (cancel_but, apply_but, ok_but):
            self.pack_start(w, 0, 0)
            w.show()
        #self.show_all()
        
def oops(message, message2=None):
    dialog = gtk.MessageDialog(flags=gtk.DIALOG_MODAL,
                               type=gtk.MESSAGE_WARNING,
                               buttons=gtk.BUTTONS_CLOSE,
                               message_format=message)
    try:
        dialog.format_secondary_text(message2)
    except AttributeError:
        print >>sys.stderr, message
        print >>sys.stderr, message2
    dialog.connect('response', lambda x, y, dialog=dialog: dialog.destroy())
    dialog.show()

class AseGuiCancelException(Exception):
    pass
        
