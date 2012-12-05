import sys

msg = "\nThe gtk python module is missing or not installed properly.\n"
msg += "Is the PYTHONPATH environment variable set correctly?\n"
msg += "Please verify your installation by running on the command line:\n"
msg += "python -c 'import gtk'\n"
msg += "\n"
msg += "This module is optional and required in order to use "
msg += "ASE's simple GUI (ag).\n"
msg += "If you don't wish to use ag ignore this error, otherwise\n"
msg += "please install the package using "
msg += "your distribution package manager, i.e.:\n"
msg += "\n"
msg += "  Debian/Ubuntu: sudo apt-get python-gtk2\n"
msg += "\n"
msg += "  OpenSUSE: yast -i python-gtk\n"
msg += "\n"
msg += "  Red Hat/Fedora: yum install pygtk2\n"
msg += "\n"
msg += "or perform manual installation, preferably as non-root user,\n"
msg += "following http://www.pygtk.org/downloads.html.\n"

if locals().get('display'):
    try:
        import gtk
    except ImportError:
        print >> sys.stderr, msg
        raise
